import pysam
import sys
import csv

"""
This script defines classes to handle a VCF file.

VCFFileProcessor processes a VCF file and creates a dictionary.
In the dictionary at each position we get a list of VariantRecord objects. Each describing a variant recorded in that
position.

- populate_dictionary() is called to populate the dictionary
- get_variant_dictionary() is called to get the dictionary
"""

# we are treating SNPs and DELs as single
SNP = 0
DEL = 0
IN = 1
HOM = 0
HET = 1
HOM_ALT = 2

VCF_QUALITY_THRESHOLD = 60

class VCFRecord:
    """
    VCFRecord class handles each record read from the VCF file.
    - Used for handling record attributes as pysam's variantRecords are not straight forward
    - Not used to be set as dictionary values as record attributes.
    """
    def __init__(self, rec):
        """
        Initialize a record object and set it's attributes.
        :param rec: A VCF record
        """
        self.rec_pos = rec.pos
        self.rec_end = rec.stop
        self.rec_qual = rec.qual if rec.qual else 10
        self.rec_genotype = self._get_record_genotype(rec)
        # FILTER field of VCF file
        self.rec_filter = list(rec.filter)[0]

        # Decides if record hom_ref or not
        self.rec_not_hom = not self.is_hom_ref(self.rec_genotype)

        self.rec_chrom = rec.chrom
        self.rec_alleles = rec.alleles
        self.rec_alts = rec.alts if rec.alts else '.'
        self.rec_ref = rec.ref

    @staticmethod
    def _get_record_genotype(record):
        """
        Get the genotype of a record. It can be fetched with s['GT'] field from pysam's API.
        :param record: A VCF record.
        :return:
        """
        gts = [s['GT'] for s in record.samples.values()]
        return gts[0]

    @staticmethod
    def is_hom_ref(genotype):
        """
        Get type of a genotype (0/0, 0/1 or 1/0)
        :param genotype: Type of GT
        :return:
        """
        is_hom = True
        gt_set = set(genotype)

        if len(gt_set) > 1:
            is_hom = False
        elif len(gt_set) == 1:
            if (0 not in gt_set) and (None not in gt_set):
                is_hom = False
        else:
            raise ValueError("INVALID GENOTYPE ENCOUNTERED" + genotype)

        return is_hom

    def __str__(self):
        """
        Print a record
        :return: String to print
        """
        return_str = str(self.rec_pos) + '\t' + str(int(self.rec_qual)) + '\t' + str(self.rec_genotype) + '\t' \
                     + str(self.rec_filter) + '\t' + str(self.rec_not_hom) + '\t' + str(self.rec_alleles) + '\t' \
                     + str(self.rec_alts)
        return return_str


class VCFFileProcessor:
    """
    Processes a file and generates a dictionary
    """
    def __init__(self, file_path):
        """
        Initialize by setting a VCF file path and initiating a dictionary
        :param file_path:
        """
        self.file_path = file_path
        self.vcf_records = None
        self.genotype_dictionary = {}
        self.total_hom = 0
        self.total_het = 0
        self.total_hom_alt = 0
        self.vcf_offset = -1

    def __str__(self):
        """
        Create a string to print the dictionary generated
        :return: String for print purpose
        """
        ret_str = "POS\tQUAL\tREF\tALT\tTYPE\tCLASS\t\n"
        for pos in self.genotype_dictionary.keys():
            for rec in self.genotype_dictionary[pos]:
                ret_str += str(rec) + "\n"

        return ret_str

    def save_positional_vcf_as_bed(self, chromosome_name, output_path_name):
        # class_code = ["SNP","IN"]
        # gt_code = ["HOM","HET","HOM_ALT"]

        with open(output_path_name, 'w') as file:
            writer = csv.writer(file, delimiter='\t')

            for pos in sorted(self.genotype_dictionary):
                # print(self.genotype_dictionary[pos])
                for variant_class in [SNP,IN]:      # Assuming that SNP and DEL have een merged into a single class
                    for variant in self.genotype_dictionary[pos][variant_class]:
                        ref, alt, gt = variant
                        row = [chromosome_name, pos, pos+1, ref, alt, gt, variant_class]
                        # print(row)
                        writer.writerow(row)

    @staticmethod
    def read_positional_vcf_from_bed(bed_path_name):
        # class_code = {"SNP":0,"IN":1}
        # gt_code = {"HOM":0,"HET":1,"HOM_ALT":2}

        positional_vcf = dict()

        with open(bed_path_name, 'r') as file:
            csv_object = csv.reader(file, delimiter='\t')

            for entry in csv_object:
                # print(entry)

                chromosome_name = entry[0]
                start = int(entry[1])
                stop = int(entry[2])
                ref = entry[3]
                alt = entry[4]
                genotype = int(entry[5])
                variant_class = int(entry[6])
                record = (ref, alt, genotype)

                # print(start, record)

                if start not in positional_vcf:
                    positional_vcf[int(start)] = [[],[],[]]

                positional_vcf[start][int(variant_class)].append(record)

        return positional_vcf

    @staticmethod
    def get_genotype_class(rec, ref, alt):
        """
        Get genotype class of a record
        :param rec: VCF record
        :param ref: Reference base
        :param alt: Alternate allele
        :return:
        """
        if len(ref) == 1 and len(alt) == 1:
            return 'SNP'
        elif len(ref) < len(alt):
            return 'IN'
        elif len(ref) > len(alt):
            return 'DEL'
        else:
            raise ValueError('INVALID GENOTYPE CLASS \n' + rec)

    @staticmethod
    def get_genotype_type(genotype):
        """
        Get type of a genotype (0/0, 0/1 or 1/0)
        :param genotype: Type of GT
        :return:
        """
        g_type = ""
        if len(genotype) == 2:
            if genotype[0] == genotype[1] and genotype[0] == 0:
                g_type = "Hom"
            elif genotype[0] != genotype[1]:
                g_type = "Het"
            else:
                g_type = "Hom_alt"
        elif len(genotype) == 1:
            if genotype[0] == 0:
                g_type = "Hom"
            if genotype[0] == 1:
                g_type = "Hom_alt"
        else:
            raise ValueError("INVALID GENOTYPE ENCOUNTERED" + genotype)
        return g_type

    def _initialize_dictionary(self, ref_pos):
        """
        Initialize a key position of the VCF dictionary
        :param ref_pos: Variant position
        :return:
        """
        if ref_pos not in self.genotype_dictionary.keys():
            self.genotype_dictionary[ref_pos] = [[], [], []]

    def _genotype_indexer(self, genotype):
        if genotype == 'Hom':
            return HOM
        if genotype == 'Het':
            return HET
        if genotype == 'Hom_alt':
            return HOM_ALT

    def _handle_delete(self, pos, ref, alt, genotype):
        """
        Process a record that has deletes.
        Deletes are usually grouped together, so we break each of the deletes to make a list.
        :param rec: VCF record containing a delete
        :return: A list of delete attributes
        """
        delete_list = []
        for i in range(0, len(ref)):
            if i < len(alt):
                continue
            ref_seq = ref[i]
            alt_seq = '*'
            pos_del = pos + i + self.vcf_offset
            delete_list.append((pos_del, ref_seq, alt_seq, genotype))
        return delete_list

    def _update_dictionary(self, pos, variant_record, genotype_class):
        """
        Add a record to the dictionary
        :param variant_record: Record to be added in the dictionary
        :return:
        """
        self._initialize_dictionary(pos)
        self.genotype_dictionary[pos][genotype_class].append(variant_record)

    def _process_genotype_by_class(self, rec, genotype_class, genotype_type, alt):
        """
        Process a genotype by class. Mostly handles different classes of genotypes differently.
        :param rec: VCF record object (VCFRecord)
        :param genotype_class: Class of the genotype
        :param genotype_type: Type of the GT
        :param alt: Alternate allele
        :return:
        """
        if genotype_class == 'SNP':
            position_based_vcf = (rec.rec_ref, alt, self._genotype_indexer(genotype_type))
            self._update_dictionary(rec.rec_pos + self.vcf_offset, position_based_vcf, SNP)
        elif genotype_class == 'DEL':
            list_of_records = self._handle_delete(rec.rec_pos, rec.rec_ref, alt, genotype_type)
            for record in list_of_records:
                pos, ref_seq, alt_seq, genotype = record
                position_based_vcf = (ref_seq, alt_seq, self._genotype_indexer(genotype))
                self._update_dictionary(pos, position_based_vcf, DEL)
        elif genotype_class == 'IN':
            if len(rec.rec_ref) > 1:
                rec.rec_ref, alt = self._trim_insert_sequences(rec.rec_ref, alt)

            position_based_vcf = (rec.rec_ref, alt, self._genotype_indexer(genotype_type))
            self._update_dictionary(rec.rec_pos + self.vcf_offset, position_based_vcf, IN)

    def _trim_insert_sequences(self, ref_seq, alt_seq):
        length = len(alt_seq) - len(ref_seq)
        trimmed_ref_seq = ref_seq[0]
        trimmed_alt_seq = alt_seq[:length+1]

        return trimmed_ref_seq, trimmed_alt_seq

    def _update_count(self, genotype_type):
        """
        Keep track of total number of variations found
        :param genotype_type: Type of genotype ('Hom', 'Het', 'Hom_alt')
        :return:
        """
        if genotype_type == 'Hom':
            self.total_hom += 1
        elif genotype_type == 'Het':
            self.total_het += 1
        else:
            self.total_hom_alt += 1

    def _parse_through_records(self, vcf_record):
        """
        Handle multiple alleles of a single record.
        :param vcf_record: VCFRecord object
        :return:
        """
        # For each allele in the record
        for alt in vcf_record.rec_alts:
            # get genotype class and type
            genotype_class = self.get_genotype_class(vcf_record, vcf_record.rec_ref, alt)
            genotype_type = self.get_genotype_type(vcf_record.rec_genotype)
            # update genotype count
            self._update_count(genotype_type)
            # process and add the genotype to the dictionary
            self._process_genotype_by_class(vcf_record, genotype_class, genotype_type, alt)

    def _get_filtered_records(self, hom_filter):
        """
        Filter records to be added to the dictionary
        :return: Filtered record list
        """
        filtered_records = []
        # For each record
        for record in self.vcf_records:
            # Create VCFRecord object for that record
            vcf_record = VCFRecord(record)

            # Filter homozygous records if hom_filter is true
            if hom_filter is True and vcf_record.rec_not_hom is False:
                continue

            # If the record can be added to the dictionary add it to the list
            if vcf_record.rec_filter == 'PASS':
                filtered_records.append(vcf_record)
            elif vcf_record.rec_qual > VCF_QUALITY_THRESHOLD:
                filtered_records.append(vcf_record)

            # print(vcf_record)

        # Return the list
        return filtered_records

    def _generate_dictionary_from_records(self, records):
        """
        Add records to the dictionary
        :param records: List of VCFRecord(s)
        :return:
        """

        # For each record
        for record in records:
            # Parse through the alleles and add it to dictionary if passes the filters
            self._parse_through_records(record)

    def get_genotype_counts(self):
        """
        Return total variation count of each genotype
        :return:
        """
        return self.total_hom, self.total_het, self.total_hom_alt

    def get_variant_dictionary(self):
        """
        Return the variant dictionary
        :return:
        """
        return self.genotype_dictionary

    def populate_dictionary(self, contig, start_pos, end_pos, hom_filter):
        """
        Process a file of a conting and site.
        :param contig: Contig (ex: chr3)
        :param site: Site (:100000-200000)
        :return:
        """

        try:
            self.vcf_records = pysam.VariantFile(self.file_path).fetch(contig, start_pos, end_pos)
        except IOError:
            sys.stderr.write("VCF FILE READ ERROR")

        # Filter the records
        filtered_records = self._get_filtered_records(hom_filter)
        # Generate dictionary
        self._generate_dictionary_from_records(filtered_records)

