"""
write metadata and output sequences at the same time
"""

from Bio import SeqIO
from datetime import datetime
from epiweeks import Week, Year
import sys
import json
import argparse
import warnings
import re
import pycountry

from datafunk.travel_history import *
from datafunk.travel_history import cities_dict, countries_list, subdivisions_dict, others
from datafunk.add_epi_week import date_string_to_epi_week

"""Don't edit these two lists please:
"""
_fields_gisaid = ['covv_accession_id', 'covv_virus_name', 'covv_location', 'covv_collection_date',
                 'covv_add_host_info', 'covv_assembly_method', 'covv_gender', 'covv_host',
                 'covv_passage', 'covv_patient_age', 'covv_seq_technology',
                 'covv_specimen', 'covv_subm_date']

_fields_edin = ['edin_header', 'edin_admin_0', 'edin_admin_1', 'edin_admin_2',
                'edin_travel', 'edin_date_stamp', 'edin_omitted', 'edin_epi_week',
                'edin_flag']

"""You can edit this list:
"""
fields = []


"""Functions
"""
def fix_gisaid_json_dict(gisaid_json_dict):
    """
    Remove commas from fields inside json dict
    """

    newDict = {}
    for x,y in gisaid_json_dict.items():
        newDict[x] = str(y).replace(',', '')

    return(newDict)


def get_admin_levels_from_json_dict(gisaid_json_dict, warnings = True):
    """
    get location strings from the gisaid location field
    use pycountry /
    """
    location_strings = [x.strip() for x in gisaid_json_dict['covv_location'].split("/")]

    while len(location_strings) < 4:
        location_strings.append("")

    continent = location_strings[0]
    country = location_strings[1]
    subdivision = location_strings[2]
    subsubdivision = location_strings[3]

    if any([x == country for x in ['England', 'Northern Ireland', 'Scotland', 'Wales']]):
        country = 'United Kingdom'
        subdivision = location_strings[1]
        subsubdivision = location_strings[2]

    if any([x == country for x in ['Alaska']]):
        country = 'USA'
        subdivision = location_strings[1]
        subsubdivision = location_strings[2]

    # some check here that there's a match to a real country
    # using pycountries?
    # First, these countries are known exceptions:
    if warnings:
        if all([country != x for x in ['Iran', 'South Korea', 'Russia', 'Korea', 'Democratic Republic of the Congo']]):
            try:
                pycountry.countries.lookup(country)
            except LookupError:
                warnings.warn('Check country flagged for ' + gisaid_json_dict['covv_accession_id'] + \
                              '  ("' + country + '")')

                if len(gisaid_json_dict['edin_flag']) == 0:
                    gisaid_json_dict['edin_flag'] = 'check_country'
                elif len(gisaid_json_dict['edin_flag']) > 0:
                    gisaid_json_dict['edin_flag'] = gisaid_json_dict['edin_flag'] + ':check_country'

    if country == 'United Kingdom':
        country = 'UK'

    if country == 'Korea':
        country = 'South Korea'

    if country == 'Democratic Republic of the Congo':
        country = 'DRC'

    gisaid_json_dict['edin_admin_0'] = country
    gisaid_json_dict['edin_admin_1'] = subdivision
    gisaid_json_dict['edin_admin_2'] = subsubdivision

    return(gisaid_json_dict)


def expand_dict(dict, fields_list_required, fields_list_optional):
    """
    check fields in a dict
    make/delete them as appropriate
    """

    for x in fields_list_required:
        if x not in dict:
            dict[x] = ''

    for x in fields_list_optional:
        if x not in dict:
            dict[x] = ''

    return(dict)


def get_csv_order_and_record_dict(csv_file, fields_list_required, fields_list_optional):
    """
    Read all info in a csv metadata file into memory.

    old_records is a nested dict with EPI_IDs as the
    top-level keys and a dict of key: value pairs
    as the top-level values
    """
    first = True
    old_records = {}
    record_order = []
    with open(csv_file, 'r') as f:
        for line in f:
            l = line.strip().split(',')
            if first:
                keys = l
                first = False
                continue

            if len(l) != len(keys):
                sys.exit('Badly formatted csv file: ' + csv_file.csv + ', exiting program')

            d = {x: y for x,y in zip(keys, l)}
            d = expand_dict(dict = d, fields_list_required = fields_list_required, fields_list_optional = fields_list_optional)

            ID = d['covv_accession_id']
            record_order.append(ID)
            old_records[ID] = d

    return(record_order, old_records)


def get_json_order_and_record_dict(json_file, fields_list_required, fields_list_optional):
    """
    Read all info in a GISAID json dump into memory.

    all_records is a nested dict with EPI_IDs as the
    top-level keys and a dict of key: value pairs
    as the top-level values
    """
    all_records = {}
    record_order = []
    with open(json_file, 'r') as f:
        for jsonObj in f:

            d = fix_gisaid_json_dict(json.loads(jsonObj))
            d = fix_seq_in_gisaid_json_dict(d)
            d = expand_dict(dict = d, fields_list_required = fields_list_required, fields_list_optional = fields_list_optional)

            ID = d['covv_accession_id']
            record_order.append(ID)
            all_records[ID] = d

    return(record_order, all_records)


def check_gisaid_date(dict):
    """

    """
    date = dict['covv_collection_date']
    regex = re.compile('\d{4}-\d{2}-\d{2}')
    match = re.search(regex, date)
    if not match:
        dict['edin_omitted'] = 'True'

        if len(dict['edin_flag']) == 0:
            dict['edin_flag'] = 'omitted_date'
        elif len(dict['edin_flag']) > 0:
            dict['edin_flag'] = dict['edin_flag'] + ':omitted_date'

    return(dict)


def check_edin_omitted_file(dict, omit_set):
    """
    update omission field to True if record
    is in omissions file
    """
    if omit_set:
        if dict['covv_accession_id'] in omit_set:
            dict['edin_omitted'] = 'True'

            if len(dict['edin_flag']) == 0:
                dict['edin_flag'] = 'omitted_file'
            elif len(dict['edin_flag']) > 0:
                dict['edin_flag'] = dict['edin_flag'] + ':omitted_file'

    return(dict)


def check_UK_sequence(gisaid_json_dict, exclude_uk = False):
    """
    Exclude UK sequences
    """
    if exclude_uk:
        header = gisaid_json_dict['edin_header']
        for country in ['England/', 'Scotland/', 'Wales/', 'Northern Ireland/']:
            if country.lower() in header.lower():
                gisaid_json_dict['edin_omitted'] = 'True'

                if len(gisaid_json_dict['edin_flag']) == 0:
                    gisaid_json_dict['edin_flag'] = 'omitted_UK'
                elif len(gisaid_json_dict['edin_flag']) > 0:
                    gisaid_json_dict['edin_flag'] = gisaid_json_dict['edin_flag'] + ':omitted_UK'

    return(gisaid_json_dict)


def update_edin_date_stamp_field(gisaid_json_dict):
    """
    date stamp the new records to write
    """
    mydate = str(datetime.date(datetime.now()))
    gisaid_json_dict['edin_date_stamp'] = mydate
    return(gisaid_json_dict)


def update_edin_epi_week_field(gisaid_json_dict):
    """
    record epi week by parsing sample collection date
    NB this will break in 2021!
    """
    if 'omitted_date' in gisaid_json_dict['edin_flag']:
        return(gisaid_json_dict)

    collection_date = gisaid_json_dict['covv_collection_date']

    last_2019 = Week(2019, 52)
    weeks = list(Year(2020).iterweeks())
    weeks.append(last_2019)

    # Rachel's function returns None if nothing found
    epi_week = date_string_to_epi_week(collection_date, weeks)

    if epi_week:
        gisaid_json_dict['edin_epi_week'] = epi_week

    return(gisaid_json_dict)


def get_one_metadata_line(dict, fields_list, sep = ','):
    """
    return a string of dict.values() formatted
    for writing in the order in fields_list

    l is a 'sep'-separated string formatted for printing
    """
    l = sep.join([dict[x] for x in fields_list]) + '\n'
    return(l)


def fix_seq_in_gisaid_json_dict(gisaid_json_dict):
    """
    strip whitespace and newline characters from the
    ['sequence'] field of a gisaid json object
    """

    def fix_seq(seq):
        newseq = re.sub(r"\s+", '', seq).replace('\n','')
        return(newseq)

    newDict = {}
    for x,y in gisaid_json_dict.items():
        if x == 'sequence':
            newDict['sequence'] = fix_seq(y)
        else:
            newDict[x] = y

    return(newDict)


def fix_header(header):
    """
    parse fasta header and remove problems
    """
    fixed_header = header.replace(' ', '_')\
        .replace("hCoV-19/","")\
        .replace("hCov-19/","")\
        .replace("PENDING", "")\
        .replace("UPLOADED", "")\
        .replace("None", "")

    return(fixed_header)


def add_header_to_json_dict(gisaid_json_dict):
    """
    make a sequence identifier from the gisaid dump json-format data.
    Function input (gisaid_json_dict) is one record from the dump
    """

    gisaid_json_dict = get_admin_levels_from_json_dict(gisaid_json_dict, warnings = False)

    myTempHeader = gisaid_json_dict['covv_virus_name'] + '|' + \
                    gisaid_json_dict['covv_accession_id'] + '||' + \
                    gisaid_json_dict['edin_admin_0'] + '|' + \
                    gisaid_json_dict['edin_admin_1'] + '|' + \
                    gisaid_json_dict['edin_admin_2'] + '|' + \
                    gisaid_json_dict['covv_collection_date']

    gisaid_json_dict['edin_header'] = fix_header(myTempHeader)

    return(gisaid_json_dict)


def parse_omissions_file(file):
    """
    Parse a file of records to omit.

    Relies on there being a regex match to "EPI_ISL_\d{6}" in a line.

    Only returns the first match to the regex
    """
    # TO DO: allow regexes to e.g. 'bat', 'pangolin' in the omissions file

    IDs = []
    regex = re.compile('EPI_ISL_\d{6}')
    file_is_fasta = file.split('.')[-1][0:2] == 'fa'

    with open(file, 'r') as f:
        for line in f:
            if file_is_fasta:
                if line[0] != '>':
                    continue
            elif line.startswith("#"):
                continue
            match = re.search(regex, line)
            if match:
                ID = match.group()
                IDs.append(ID)

    return(IDs)


def write_metadata_output(output,
                  new_records_list,
                  new_records_dict,
                  old_records_list,
                  old_records_dict,
                  fields_list):

    """
    write a csv-format outfile to file
    """
    out = open(output, 'w')

    out.write(','.join(fields_list) + '\n')

    for record in old_records_list:
        do = old_records_dict[record]
        lo = get_one_metadata_line(dict=do, fields_list=fields_list)
        out.write(lo)

    for record in new_records_list:
        dn = new_records_dict[record]
        ln = get_one_metadata_line(dict=dn, fields_list=fields_list)
        out.write(ln)

    out.close()
    pass


def write_fasta_output(output,
                  new_records_list,
                  new_records_dict,
                  old_records_list,
                  old_records_dict):
    if output:
        out = open(output, 'w')
    else:
        out = sys.stdout

    for record in old_records_list:
        if old_records_dict[record]['edin_omitted'] == 'True':
            continue
        else:
            header = old_records_dict[record]['edin_header']
            seq = old_records_dict[record]['sequence']
            out.write('>' + header + '\n')
            out.write(seq + '\n')

    for record in new_records_list:
        if new_records_dict[record]['edin_omitted'] == 'True':
            continue
        else:
            header = new_records_dict[record]['edin_header']
            seq = new_records_dict[record]['sequence']
            out.write('>' + header + '\n')
            out.write(seq + '\n')

    if output:
        out.close()
    pass


"""Program
"""
def process_gisaid_data(input_json,
                        input_omit_file_list,
                        input_metadata,
                        output_fasta,
                        output_metadata,
                        exclude_uk=False,
                        exclude_undated=False):

    # logfile = open(output + '.log', 'w')
    if input_omit_file_list:
        temp = []
        for file in input_omit_file_list:
            temp.extend(parse_omissions_file(file))

        omitted_IDs = set(temp)
    else:
        omitted_IDs = False


    if input_metadata != 'False':
        # Check that all required fields were in the csv file
        csv_header = next(open(input_metadata, 'r')).strip().split(',')
        if not all([x in csv_header for x in _fields_gisaid + _fields_edin]):
            sys.exit('There were missing mandatory fields in ' + input_metadata)

        # add optional fields from the old csv file to the output
        for x in csv_header:
            if x not in _fields_gisaid + _fields_edin:
                if x not in fields:
                    fields.append(x)


        old_records = get_csv_order_and_record_dict(input_metadata,
                                                    fields_list_required = _fields_gisaid + _fields_edin,
                                                    fields_list_optional = fields)

        old_records_list = old_records[0]
        old_records_dict = old_records[1]

    else:
        old_records_list = []
        old_records_dict = {}


    all_records = get_json_order_and_record_dict(input_json, \
                                                fields_list_required = _fields_gisaid + _fields_edin, \
                                                fields_list_optional = fields)

    all_records_list = all_records[0]
    all_records_dict = all_records[1]

    new_records_list = [x for x in all_records_list if x not in set(old_records_list)]
    new_records_dict = {x: all_records_dict[x] for x in all_records_list if x not in set(old_records_list)}

    # update date stamp field
    new_records_dict = {x: update_edin_date_stamp_field(all_records_dict[x]) for x in new_records_dict.keys()}

    # update admin level
    new_records_dict = {x: get_admin_levels_from_json_dict(all_records_dict[x]) for x in new_records_dict.keys()}

    # include a header field in each dictionary (just for writing the fasta file):
    new_records_dict = {x: add_header_to_json_dict(all_records_dict[x]) for x in new_records_dict.keys()}

    # get travel history
    new_records_dict = {x: get_travel_history(all_records_dict[x]) for x in new_records_dict.keys()}

    # check gisaid collection date formatted correctly
    new_records_dict = {x: check_gisaid_date(all_records_dict[x]) for x in new_records_dict.keys()}

    # if gisaid collection date formatted correctly, we can add epi week
    new_records_dict = {x: update_edin_epi_week_field(all_records_dict[x]) for x in new_records_dict.keys()}

    # check if sequence is in omissions file
    new_records_dict = {x: check_edin_omitted_file(all_records_dict[x], omitted_IDs) for x in new_records_dict.keys()}

    if exclude_uk:
        # check if sequence is from the UK
        new_records_dict = {x: check_UK_sequence(all_records_dict[x], exclude_uk = True) for x in new_records_dict.keys()}


    if output_metadata:
        write_metadata_output(output = output_metadata,
                              new_records_list = new_records_list,
                              new_records_dict = new_records_dict,
                              old_records_list = old_records_list,
                              old_records_dict = old_records_dict,
                              fields_list = _fields_edin + fields + _fields_gisaid)


    write_fasta_output(output = output_fasta,
                       new_records_list = new_records_list,
                       new_records_dict = new_records_dict,
                       old_records_list = old_records_list,
                       old_records_dict = old_records_dict)



    # logfile.close()
    pass




















#