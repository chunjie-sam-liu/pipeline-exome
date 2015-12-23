'''
Created on 08/01/2014

@author: antonior
'''

class region():
    '''
    This class manages chromosomic regions
    '''


    def __init__(self,chrom,start,end,value=None):
        '''
        Constructor
        '''
        self.chrom=chrom
        self.start=int(start)
        self.end=int(end)
        if value!=None:
            self.value=int(value)
        else:
            self.value=value
    def __str__ (self):
        
        if self.value!=None:
            return self.chrom+' '+str(self.start)+' '+str(self.end)+' '+str(self.value)+'\n' 
    
        else:
            return self.chrom+' '+str(self.start)+' '+str(self.end)+'\n'
        
    def overlap_type (self,other):
        
        if self.chrom!=other.chrom:
            return 0
        
        elif other.start>=self.end:
            return -1
        
        elif other.end<=self.start:
            return 0
        
        elif other.start<=self.end and other.start>self.start and other.end>self.end:
            return 1
        
        elif other.start<=self.start and other.end>self.start and other.end<self.end:
            return 2
        
        elif other.start>self.start and other.end<self.end:
            return 3
        
        elif other.start<=self.start and other.end>=self.end:
            return 4
        
        else:
            raise ValueError()


    def __sub__ (self,other):
            
        overlap=self.overlap_type(other)
            
        if overlap==0 or overlap==-1:
            r1=self
            r=[r1]
            return r
        if overlap==1:
            r1=region(self.chrom,self.start,other.start,self.value)
            r=[r1]
            return r
        if overlap==2:
            r1=region(self.chrom,other.end,self.end,self.value)
            r=[r1]
            return r
        if overlap==3:
            r1=region(self.chrom,self.start,other.start,self.value)
            r2=region(self.chrom,other.end,self.end,self.value)
            r=[r1,r2]
            return r
        if overlap==4:
            r1=[]
            return r1
                
         
        