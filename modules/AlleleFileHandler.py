import json

class AlleleFileHandler:
    '''
    Opens a JSON file containing a dictionary with key = sites and value = list of observed alleles, each with relevant
    attributes, E.g.:
                {
                    "allele_sequence": "TGG",
                    "base_qualities": [34,34,34],
                    "map_quality": 60,
                    "read_direction": true,
                    "read_id": "HSQ1004:134:C0D8DACXX:2:2103:7674:57255"
                },

    and loads the data into a dictionary.
    '''

    def __init__(self, allele_file_path):
        self.file_path = allele_file_path
        self.allele_dictionary = None

        self.load_from_file()

    def load_from_file(self):
        '''
        Load the JSON file into a dictionary
        '''
        with open(self.file_path, 'r') as file:
            self.allele_dictionary = json.load(file)

    def get_allele_site_dictionary(self):
        '''
        Get the allele_dictionary containing all sites, with full list of alleles sampled at each site
        :return: allele_dictionary
        '''
        return self.allele_dictionary
