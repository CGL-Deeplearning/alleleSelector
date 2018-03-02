DEFAULT_MIN_MAP_QUALITY = 5
IMAGE_HEIGHT = 300
IMAGE_WIDTH = 300
IMAGE_BUFFER = 0
CIGAR_MATCH = 0
CIGAR_IN = 1
CIGAR_DEL = 2
MAX_COLOR_VALUE = 254.0
BASE_QUALITY_CAP = 40.0
MAP_QUALITY_CAP = 60.0
MAP_QUALITY_FILTER = 10.0


class ImageChannels:
    """
    Handles how many channels to create for each base and their way of construction.
    """
    @staticmethod
    def get_base_color(base):
        """
        Get color based on a base.
        - Uses different band of the same channel.
        :param base:
        :return:
        """
        if base == 'A':
            return 250.0
        if base == 'C':
            return 100.0
        if base == 'G':
            return 180.0
        if base == 'T':
            return 30.0
        if base == '*' or 'N':
            return 5.0

    @staticmethod
    def get_base_quality_color(base_quality):
        """
        Get a color spectrum given base quality
        :param base_quality: value of base quality
        :return:
        """
        c_q = min(base_quality, BASE_QUALITY_CAP)
        color = MAX_COLOR_VALUE * c_q / BASE_QUALITY_CAP
        return color

    @staticmethod
    def get_map_quality_color(map_quality):
        """
        Get a color spectrum given mapping quality
        :param map_quality: value of mapping quality
        :return:
        """
        c_q = min(map_quality, MAP_QUALITY_CAP)
        color = MAX_COLOR_VALUE * c_q / MAP_QUALITY_CAP
        return color

    @staticmethod
    def get_strand_color(is_rev):
        """
        Get color for forward and reverse reads
        :param is_rev: True if read is reversed
        :return:
        """
        if is_rev is True:
            return 240
        else:
            return 70

    @staticmethod
    def get_match_ref_color(is_match):
        """
        Get color for base matching to reference
        :param is_match: If true, base matches to reference
        :return:
        """
        if is_match is True:
            return MAX_COLOR_VALUE * 0.2
        else:
            return MAX_COLOR_VALUE * 1.0

    @staticmethod
    def get_alt_support_color(is_in_support):
        """
        ***NOT USED YET***
        :param is_in_support:
        :return:
        """
        if is_in_support is True:
            return MAX_COLOR_VALUE * 1.0
        else:
            return MAX_COLOR_VALUE * 0.6

    @staticmethod
    def get_empty_channels():
        """
        Get empty channel values
        :return:
        """
        return [0, 0, 0, 0, 0, 0]

    @staticmethod
    def get_channels(attribute_tuple):
        """
        Get a bases's channel construction
        :return: [color spectrum of channels based on base attributes]
        """
        base, base_q, map_q, is_rev, is_match, is_supporting = attribute_tuple
        base_color = ImageChannels.get_base_color(base)
        base_quality_color = ImageChannels.get_base_quality_color(base_q)
        map_quality_color = ImageChannels.get_map_quality_color(map_q)
        strand_color = ImageChannels.get_strand_color(is_rev)
        match_color = ImageChannels.get_match_ref_color(is_match)
        get_support_color = ImageChannels.get_alt_support_color(is_supporting)

        return [base_color, base_quality_color, map_quality_color, strand_color, match_color, get_support_color]

    @staticmethod
    def get_ref_channels(base):
        """
        Get a reference bases's channel construction
        :param base: Reference base
        :return: [color spectrum of channels based on some default values]
        """
        base_color = ImageChannels.get_base_color(base)
        base_quality_color = ImageChannels.get_base_quality_color(60)
        map_quality_color = ImageChannels.get_map_quality_color(60)
        strand_color = ImageChannels.get_strand_color(is_rev=False)
        match_color = ImageChannels.get_match_ref_color(is_match=True)
        support_color = ImageChannels.get_alt_support_color(is_in_support=True)

        return [base_color, base_quality_color, map_quality_color, strand_color, match_color, support_color]

    # RGB image creator
    # ---ONLY USED FOR TESTING--- #
    @staticmethod
    def get_empty_rgb_channels():
        return [0, 0, 0, 255]

    @staticmethod
    def get_color_for_base_rgb(ref, base):
        if ref == base and ref != '*':
            return 255, 255, 255
        elif base == 'A':
            return 255, 0, 0
        elif base == 'C':
            return 255, 255, 0
        elif base == 'T':
            return 0, 0, 255
        elif base == 'G':
            return 0, 255, 0
        else:
            return 255, 0, 255

    @staticmethod
    def get_channels_only_rgb(attribute_tuple, ref_base):
        base, base_q, map_q, is_rev, is_match, is_supporting = attribute_tuple
        base_color = ImageChannels.get_base_color(base)
        base_quality_color = ImageChannels.get_base_quality_color(base_q)
        map_quality_color = ImageChannels.get_map_quality_color(map_q)
        strand_color = ImageChannels.get_strand_color(is_rev)
        match_color = ImageChannels.get_match_ref_color(is_match)
        support_color = ImageChannels.get_alt_support_color(is_supporting)
        r, g, b = ImageChannels.get_color_for_base_rgb(ref_base, base)

        return [r, g, b, support_color]

    @staticmethod
    def get_ref_channels_rgb(base):
        base_color = ImageChannels.get_base_color(base)
        base_quality_color = ImageChannels.get_base_quality_color(60)
        map_quality_color = ImageChannels.get_map_quality_color(60)
        strand_color = ImageChannels.get_strand_color(is_rev=False)
        get_match_color = ImageChannels.get_match_ref_color(is_match=True)
        r, g, b = ImageChannels.get_color_for_base_rgb('', base)
        support_color = ImageChannels.get_alt_support_color(is_in_support=True)

        return [r, g, b, support_color]