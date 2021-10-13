
from collections import namedtuple
from pyteomics import mzxml
from threading import Thread
import urllib.request
import queue

Fragment = namedtuple('Fragment', ['Weight', 'Abundance'], defaults=[0.0, 0])
Precursor = namedtuple('Precursor', ['Id', 'PrecursorWeight', 'PrecursorCharge', 'Fragments'], defaults=['', 0.0, 0, {}])
ReferenceMatch = namedtuple('ReferenceMatch',['ProteinName','StartIndex','AminoAcids'], defaults=['',0,{}])
FragmentReferenceMatches = namedtuple('FragmentReferenceMatch',['IonCharge','ChargeAmount','Fragment','ReferenceMatches'],defaults=['',1,None,{}])

def generate_precursor_from_content(content):
    id = content['id']
    precursor_charge = content['precursorMz'][0]['precursorCharge']
    precursor_weight = content['precursorMz'][0]['precursorMz']
    fragments_weights = content['m/z array']
    fragments_abundances = content['intensity array']
    raw_fragments = zip(fragments_weights,fragments_abundances)        
    fragments = []
    for raw_fragment in raw_fragments:
        fragment = Fragment(raw_fragment[0],raw_fragment[1])
        fragments.append(fragment)
    return Precursor(id,precursor_weight,precursor_charge,fragments)

def generate_precursors_from_mzml_file(file_path):
    file_contents = mzxml.read(file_path)
    contents = []
    for file_content in file_contents:
        contents.append(file_content)
    observations = []        
    for content in contents:
        observation = generate_precursor_from_content(content)
        observations.append(observation)
    return observations       

def print_precursors(precursors):
    number_of_precursors = len(precursors)
    total_number_of_fragments = 0
    lowest_number_of_fragments = 999999
    largest_number_of_fragments = 0
    for precursor in precursors:
        number_of_fragments = len(precursor.Fragments)
        if number_of_fragments < lowest_number_of_fragments: lowest_number_of_fragments = number_of_fragments
        if number_of_fragments > largest_number_of_fragments: largest_number_of_fragments = number_of_fragments
        total_number_of_fragments += number_of_fragments
    average_number_of_fragments = total_number_of_fragments/number_of_precursors
    
    formatted_number_of_precursors = str(number_of_precursors)
    formatted_total_number_of_fragments = str(total_number_of_fragments)
    formatted_average_number_of_fragments = str(average_number_of_fragments)
    formatted_lowest_number_of_fragments = str(lowest_number_of_fragments)
    formatted_largest_number_of_fragments = str(largest_number_of_fragments)

    print('number of precursors in MzXML file:' + formatted_number_of_precursors)
    print('total number of fragments in MzXML File: ' + formatted_total_number_of_fragments)    
    print('average number of fragments per precursor: ' + formatted_average_number_of_fragments)
    print('lowest number of fragments: ' + formatted_lowest_number_of_fragments) 
    print('highest number of fragments: ' + formatted_largest_number_of_fragments) 

def generate_charge_amount_adjusted_fragments(fragments, charge_amount):
    adjusted_fragments = []
    for fragment in fragments:
        adjusted_fragment = Fragment(fragment.Weight * charge_amount, fragment.Abundance)
        adjusted_fragments.append(adjusted_fragment)
    return adjusted_fragments

def generate_fragments_from_precursors(precursors,charge_amount):
    fragments = []
    for precursor in precursors:
        fragment = Fragment(precursor.PrecursorWeight * charge_amount, 999999)
        fragments.append(fragment)
        adjusted_fragments = generate_charge_amount_adjusted_fragments(precursor.Fragments,charge_amount)
        fragments.extend(adjusted_fragments)
    return fragments     

def perform_parallel_web_requests(addresses, no_workers):
    class Worker(Thread):
        def __init__(self, request_queue):
            Thread.__init__(self)
            self.queue = request_queue
            self.results = []

        def run(self):
            while True:
                content = self.queue.get()
                if content == "":
                    break
                request = urllib.request.Request(content)
                response = urllib.request.urlopen(request)
                self.results.append(response.json())
                self.queue.task_done()
    q = queue.Queue()
    for url in addresses:
        q.put(url)
    workers = []
    for _ in range(no_workers):
        worker = Worker(q)
        worker.start()
        workers.append(worker)
    for _ in workers:
        q.put("")
    for worker in workers:
        worker.join()
    r = []
    for worker in workers:
        r.extend(worker.results)
    return r

def generate_url_from_fragment(fragment, ion_charge, ppm_tolerance):
    return 'http://hypedsearchservice.azurewebsites.net/api/proteinmatch?ion_charge=' + ion_charge + '&weight=' + str(fragment.Weight) + '&ppm_tolerance=' + str(ppm_tolerance)
    
def generate_reference_matches_from_json(json):
    reference_matches = []
    for item in json:
        protein_name = item['ProteinName']
        start_index = item['StartIndex']
        kmers = item['KMers']
        reference_match = ReferenceMatch(protein_name,start_index,kmers)
        reference_matches.append(reference_match)
    return reference_matches

def generate_reference_matches_from_fragment(fragment, ion_charge, ppm_tolerance):
    uri = 'http://hypedsearchservice.azurewebsites.net/api/proteinmatch?ion_charge=' + ion_charge + '&weight=' + str(fragment.Weight) + '&ppm_tolerance=' + str(ppm_tolerance)
    response = urllib.requests.get(uri)
    response_json = response.json()
    reference_matches = []
    for item in response_json:
        protein_name = item['ProteinName']
        start_index = item['StartIndex']
        kmers = item['KMers']
        reference_match = ReferenceMatch(protein_name,start_index,kmers)
        reference_matches.append(reference_match)
    return reference_matches

def generate_fragment_reference_matchess(fragments,ion_charge, charge_amount, ppm_tolerance):
    all_fragment_reference_matches = []
    for fragment in fragments:
        reference_matches = generate_reference_matches_from_fragment(fragment,ion_charge, ppm_tolerance)
        fragment_reference_matches = FragmentReferenceMatches(ion_charge,charge_amount,fragment,reference_matches)
        all_fragment_reference_matches.append(fragment_reference_matches)
        return all_fragment_reference_matches

def parallel_generate_fragment_reference_matchess(fragments,ion_charge, charge_amount, ppm_tolerance):
    all_fragment_reference_matches = []
    for fragment in fragments:
        urls = generate_url_from_fragment(fragment,ion_charge, ppm_tolerance)
        jsons = perform_parallel_web_requests(urls,16)
        for json in jsons:
            reference_matches = generate_reference_matches_from_json(json)
            fragment_reference_matches = FragmentReferenceMatches(ion_charge,charge_amount,fragment,reference_matches)
            all_fragment_reference_matches.append(fragment_reference_matches)
        return all_fragment_reference_matches

def generate_all_fragment_reference_matches(precursors,ion_charges,charge_amounts,ppm_tolerance):
    all_fragment_reference_matches = []
    for ion_charge in ion_charges:
        for charge_amount in charge_amounts:
            fragments = generate_fragments_from_precursors(precursors,charge_amount)
            fragment_reference_matchess = parallel_generate_fragment_reference_matchess(fragments,ion_charge,charge_amount,ppm_tolerance)
            all_fragment_reference_matches.append(fragment_reference_matchess)
    return all_fragment_reference_matches

def print_all_fragment_reference_matches(fragmentReferenceMatches):
    total_number_of_fragments = len(fragmentReferenceMatches)
    formatted_total_number_of_fragments = str(total_number_of_fragments)
    print('Total number of fragment matches: ' + formatted_total_number_of_fragments)

def countby(matches):
    result = {}
    for match in matches: 
        key = match.ProteinName
        if key in result:
            result[key] += 1
        else: 
            result[key] = 1
    return result

def sortDictionary(dict):
    sorted_values = sorted(dict.values())
    sorted_dict = {} 
    for i in sorted_values:
        for k in dict.keys():
            if dict[k] == i:
                sorted_dict[k] = dict[k]
                break
    return sorted_dict

path = '/Users/jamesdixon/Documents/Hypedsearch_All/Parse MzXML File/BM_SEC_F3_RA_AspN.mzXML'
ion_charges = ['B','Y']
charge_amounts = [1,2]
ppm_tolerance = 0.01
precursors = generate_precursors_from_mzml_file(path)
print_precursors(precursors)
all_fragment_reference_matches = generate_all_fragment_reference_matches(precursors,ion_charges,charge_amounts,ppm_tolerance)



#first_matches = first.ReferenceMatches
#counted = countby(first_matches)
#sorted = sortDictionary(counted)
#print(sorted)

#this should be a unit(really integration) test :-)
#fragment = Fragment(218.07,1000)
#reference_matches = generate_reference_matches_from_fragment(fragment,'B',.01)
#fragment_reference_matches = FragmentReferenceMatches(fragment,reference_matches)



