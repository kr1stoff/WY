#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Chaobo Ren
# @Date:   2022/3/28 14:30
# @Last Modified by:   Ming
# @Last Modified time: 2022/3/28 14:30
import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

class cgMLST(object):
    """
    cgMLST Pipeline
    """
    def __init__(self):
        """
        Init the object
        """
        pass

    def filterGenome(self,f_keep):
        """
        Filter the Genome file by their name

        :param f_keep: The file contain genome files to keep
        """
        pass

    def trainModel(self,f_ref):
        """
        Use the reference genome to train a gene model

        :param f_ref: The genome file as reference to train model
        """
        pass

    def wgMLST(self):
        """
        Use the reference to get the whole genome MLST schema
        """
        pass

    def alleleCalling(self):
        """
        Allele calling
        """
        pass

    def cgMLST(self):
        """
        Get the core genome MLST schema
        """
        pass



