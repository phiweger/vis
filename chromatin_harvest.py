from __future__ import division
import HTSeq
from collections import defaultdict, namedtuple 
from itertools import chain
import numpy as np
import pandas as pd



class Feature(object):

    def __init__(self, mydict):
        for key in mydict:
            setattr(self, key, mydict[key])

    def setval(self, name, value):
        '''
        useful when setting new attributes like "has_DMR"
        '''
        setattr(self, name, value)

    def __getitem__(self, name):
        '''
        enables feature['name'] syntax
        '''
        return getattr(self, name)



class Switch(Feature):

    def __repr__(self):
        return str(self.iv)
        # Since self.iv is a GenomicInterval object from HTSeq, it falls
        # back to its own __repr__ method.



class Dmr(Feature):

    def __repr__(self):
        return '({}, signal: {})'.format(str(self.iv), self.diffsignal)
        # Since self.iv is a GenomicInterval object from HTSeq, it falls
        # back to its own __repr__ method.



class Tfbs(Feature):

    def __repr__(self):
        return str(self.iv)



class Compartment(Feature):

    def __repr__(self):
        return str(self.iv)
        # Since self.iv is a GenomicInterval object from HTSeq, it falls
        # back to its own __repr__ method.



# The sow function creates data "fields" which we'll "harvest" later
# on. Each field represents a dimension in the analysis, such as chromatin
# state transitions ("switch"), DMR, methylation calling, expression etc.
def sow_field(fp):
    '''
    This function will return a GenomicArrayOfSets object loaded with
    any GenomicFeature. The latter is equivalent to a line in a bed 
    file.

    field = sow(fp)
    '''
    gas = HTSeq.GenomicArrayOfSets('auto', stranded=False)
    field = HTSeq.BED_Reader(fp)
    for line in field:
        gas[line.iv] += line
    return gas



def sow_expression(fp):
    '''
    Expression data as obtained by DESeq2.
    '''
    names=('name baseMean log2FoldChange ' 
    'lfcSE stat pvalue padj'.split(' '))
    df = pd.read_table(fp, header=None, names=names)
            
    grp = df.groupby('name', group_keys=False)
    df = grp.filter(lambda x: len(x) == 1)
    return df



def sow_methylation(fp):
    '''
    Loads the methylation calling and saves them in an GenomicArray.
    Note that we will not save the methylation into objects of type
    GenomicPosition storing them in GenomicArrayOfSets because this 
    would add little to convenience but a lot to memory consumption.
    '''
    ga = HTSeq.GenomicArray('auto', stranded=False, typecode='d') # double
    methylation = HTSeq.BED_Reader(fp)
    for m in methylation:
        ga[m.iv] = float(m.name)
    return ga



def sow_cpg(fp):
    '''
    Loads a track of positions that are either C or G in the CpG context,
    i.e. their occuring in each other's presence. The data is saved to a 
    GenomicArray.
    Note that we will not save the methylation into objects of type
    GenomicPosition storing them in GenomicArrayOfSets because this 
    would add little to convenience but a lot to memory consumption.
    '''
    ga = HTSeq.GenomicArray('auto', stranded=False, typecode='d') # double
    cpg = HTSeq.BED_Reader(fp)
    for c in cpg:
        ga[c.iv] = int(c.name) # there are only 1s in the file
    return ga



# The workhorse of the analysis, the Harvest class "harvests" all the information
# from the various data fields. It is designed to be an interactive tool, that
# facilitates the access to different data fields. 
class Harvest(object):

    def __init__(self, gene):
        self.id, self.name = gene.name.split(';')
        self.iv = gene.iv # will hold a GenomicInterval object from HTSeq
        
        self.expression = None 
        self.switch, self.dmr, self.compartment, self.tfbs = [
            defaultdict(list) for i in range(4)]

    def __repr__(self):
        return 'name: {}, id: {}'.format(self.name, self.id)



    # beforehand we need to call the sow(fp) function
    def harvest_field(self, gas, feature, overlap=False): 
        '''
        This method defines the interface to the respective fields.
        '''

        count = 0

        if overlap == True:
            if feature == 'tfbs':
                v = gas[self.iv]
                bag = set()
                for i in v.values(): # loop through positions of ChromVector
                    if i: # i.e. if the set is not empty
                        for element in i: # loop through set entries
                            bag.add(element) # bag contains GenomicFeature
                
                count = len(bag)

                for k in bag:
                    d = self._get_feature_overlapping(k)
                    obj = Tfbs(d)
                    self.tfbs[obj.mykey].append(obj)
                
        else:
            for step in gas[self.iv].steps():
                if step[1]:
                # step[1] is of type set, evaluates to False if empty;
                # get_dmr() will throw an error if called on an empty set
                    
                    if feature == 'dmr':
                        d = self._get_dmr(step)
                        obj = Dmr(d)
                        self.dmr[obj.mykey].append(obj) # type .. hyper- or hypomethylated
                        count += 1
            
                    if feature == 'switch':
                        d = self._get_feature(step)
                        obj = Switch(d)
                        self.switch[obj.mykey].append(obj)
                        count += 1
    
                    if feature == 'compartment':
                        d = self._get_feature(step)
                        obj = Compartment(d)
                        self.compartment[obj.mykey].append(obj)
                        count += 1
                    
                    # historic
                    # if feature == 'tfbs':
                        # d = self._get_feature(step)
                        # obj = Tfbs(d)
                        # self.tfbs[obj.mykey].append(obj)
                        # count += 1

        setattr(self, 'num_' + feature , count)



    def harvest_expression(self, e): # e .. expression
        '''
        Creates a dict self.expression by looking up the genes name in 
        the lookup table generated with sow_expression().
        '''
        try:
            self.expression =  e[e.name.str.match(
                # .match() is more strict than .contains()
                '^' + self.name + '$' 
                # otherwise e.g. SAMD11 and SAMD11P1 not distinguished
                )].to_dict(orient='record')[0]

            # will return a dict with empty list if the gene is not found
            # in the file with the differential expression; the [0] is needed
            # because with the 'record' option pandas returns a list[dict], so
            # we extract the dict
        except IndexError:
            print('No expression data available for gene ' + self.name + '. '
                'Sad, is it not?')



    #                          -------------------
    # Various _get_x() methods that modify the feature in the
    # GenomicArrayOfSets according to the needs of the respective
    # information, e.g. by deciding which key is used in the 
    # defaultdict that holds the features.

    def _get_feature_overlapping(self, feature):
        '''
        Method to turn a genomic feature in an overlapping feature set into
        a feature object, e.g. of type Tfbs. Note that while 
        _get_feature() takes a step object as input, 
        _get_feature_overlapping requires a GenomicFeature.
        '''

        d = dict([
        ('mykey', feature.name),
        ('iv', feature.iv),
        # for convenience: now we can access a GenomicArrayOfSets with
        # switch = Switch(mydict)
        # gas[switch.iv]
        ('feature', feature)
        ])

        return d



    def _get_feature(self, step): # only if non-overlapping features
        '''
        Extracts some information from steps of a GenomicArrayOfSets. This
        is the method to call for the common bed format, i.e. when all we
        need as information is the 4th column holding some annotation, e.g.
        for chromatin state transitions or interaction compartments.
        '''
        iv, feature = step[0], next(iter(step[1])) 
        # alternatives:
        # name = list(name)[0]
        # [name] = name
        # (name,) = name
        # source: stackoverflow, 6091922

        # distances
        if self.iv.strand == '+':
            dist = iv.start - self.iv.start
            dist_rel = dist / self.iv.length
        elif self.iv.strand == '-':
            dist = self.iv.start_d - iv.end
            dist_rel = dist / self.iv.length
        else:
            dist = 'NA'
            dist_rel = 'NA' 

        # This would be wrong, as we need to coonsider strandedness:
        # dist = iv.start - self.iv.start
        # dist_rel = dist / self.iv.length
        
        coverage = iv.length / self.iv.length
        # from future import __division__
        # otherwise cast to float()
        
        d = dict([
        ('mykey', feature.name),
        ('iv', iv),
        # for convenience: now we can access a GenomicArrayOfSets with
        # switch = Switch(mydict)
        # gas[switch.iv]
        ('dist', dist), 
        ('dist_rel', dist_rel), 
        ('coverage', coverage),
        ('feature', feature)
        ])
        # Note: The feature also has an attribute called "iv". The difference
        # between iv and feature.iv is that the former returns the interval
        # associated with the step, i.e. truncated DMR intervals, while the 
        # latter returns the DMRs' original intervals. This implementation 
        # potentially allows to calculate the mean methylation for the whole
        # DMR instead of the truncated one.

        # In [203]: tmp.iv
        # Out[203]: <GenomicInterval object 'chr1', [150,159), strand '.'>
        # In [204]: tmp.feature.iv
        # Out[204]: <GenomicInterval object 'chr1', [150,170), strand '.'>
        return d



    def _get_dmr(self, step):
        iv, feature = step[0], next(iter(step[1])) 
        # step[0] holds a GenomicInterval object,
        # step[1] holds a set containing a GenomicFeature object (the
        # next(iter()) construction first makes the set iterable and 
        # then grabs the first item, thus "extracting" the object from the 
        # set)
        difftype = 'hyper' if float(feature.name) > 0 else 'hypo'

        # distances
        if self.iv.strand == '+':
            dist = iv.start - self.iv.start
            dist_rel = dist / self.iv.length
        elif self.iv.strand == '-':
            dist = self.iv.start_d - iv.end
            dist_rel = dist / self.iv.length
        else:
            dist = 'NA'
            dist_rel = 'NA'         

        coverage = iv.length / self.iv.length

        d = dict([
        ('mykey', difftype),
        ('iv', iv),
        # for convenience: now we can access a GenomicArrayOfSets with
        # switch = Switch(mydict)
        # gas[switch.iv]
        ('diffsignal', float(feature.name)),
        ('dist', dist), 
        ('dist_rel', dist_rel), 
        ('coverage', coverage),
        ('feature', feature)
        ])
        return d



    #                          -------------------
    # method to set Switch.has_dmr 

    def has_dmr(self, feature):
        '''
        Here we query the switch and DMR dict and annotate each switch
        as either being covered by a DMR or not. Note that the overlap
        takes into account the full length of both features, not just the
        fragment that covers the gene (or tss).
        '''      
        # if feature == 'switch':
        #     if bool(self.switch) & bool(self.dmr):
        #         for s in chain.from_iterable(self.switch.values()): 
        #         # s .. Switch object
        #             for d in chain.from_iterable(self.dmr.values()):
        #                 if s.feature.iv.overlaps(d.feature.iv):
        #                     s.setval('has_dmr', True)
        #                 else:    
        #                     s.setval('has_dmr', False)


        if hasattr(self, 'num_' + feature) & hasattr(self, 'num_dmr'):
            
            for s in chain.from_iterable(
                getattr(self, feature).values()
                ): 
            # s .. Switch object

                # in case the attribute dmr is empty
                if not getattr(self, 'dmr'):
                    s.setval('has_dmr', False)
                else:
                    for d in chain.from_iterable(self.dmr.values()):
                        if s.feature.iv.overlaps(d.feature.iv):
                            s.setval('has_dmr', True)
                        else:    
                            s.setval('has_dmr', False)
        else:
            print 'Please harvest both - feature and DMR - first.'



    #                          -------------------
    # method to set Feature.methyl_base and  Feature.methyl_diff

    def get_methylation(self, ga, type, feature): # e.g. (ga, 'base', 'switch')
        # ga .. object GenomicArray
        '''
        Will extract baseline and differential methylation rate for 
        the feature's interval (to allow calculation e.g. og the mean
        methylation rate of a switch, DMR or TFBS). Note that the methylation
        info is returned for the original interval of the feature, not
        the one truncated by e.g. the gene's boundaries. 
        '''

        for s in chain.from_iterable(
            getattr(self, feature).values()
            ):
            arr = tuple(ga[s.feature.iv])
            
            if type == 'baseline':
                s.setval('methylation_baseline', arr)
                s.setval('methylation_baseline_mean', np.mean(arr))
                s.setval('methylation_baseline_var', np.var(arr))
            elif type == 'difference':
                s.setval('methylation_difference', arr)
                s.setval('methylation_difference_mean', np.mean(arr))
                s.setval('methylation_difference_var', np.var(arr))
            else: 
                print('Permissible values for type are "baseline" and '
                    '"difference".')
                break

    #                          -------------------
    # method to set Feature.methyl_base and  Feature.methyl_diff

    def get_cpg(self, ga, feature): # e.g. (ga, 'base', 'switch')
        # ga .. object GenomicArray
        '''
        Will extract baseline and differential methylation rate for 
        the feature's interval (to allow calculation e.g. og the mean
        methylation rate of a switch, DMR or TFBS). Note that the methylation
        info is returned for the original interval of the feature, not
        the one truncated by e.g. the gene's boundaries. 
        '''

        for s in chain.from_iterable(
            getattr(self, feature).values()
            ):
            
            arr = tuple(ga[s.feature.iv])

            s.setval('cpg', arr)
            s.setval('cpg_mean', np.mean(arr))
            s.setval('cpg_var', np.var(arr))
        








    #                          -------------------
    # methylation related



    #                          -------------------
    # method to calculate a feature's methylation 



    #                          -------------------
    # differential exon usage related?



# include:
# DEXSeq? (e.g. position of DEU exon)
# compartments from arrays

