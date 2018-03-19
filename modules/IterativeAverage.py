

class IterativeAverage:
    """
    A class that keeps track of the average of all inputs, by calculating the iterative average (not storing a list)
    """
    def __init__(self):
        self.average = 0
        self.n = 0

    def update(self, x):
        """
        Update the running average with another number
        :param x:
        :return:
        """
        self.average += (1/(self.n+1))*(x-self.average)
        self.n += 1

    def get_average(self):
        """
        Return the running average
        :return: average of all previous inputs
        """
        return self.average

    def get_number_of_elements(self):
        """
        Return the total number of previous inputs
        :return: number of all previous inputs
        """
        return self.n

