# -*- coding: utf-8 -*-

"""This function cleans codon alignments by gaps and stop codons removal, allowing
	the user to select zones, generate a protein alignment and select output name and
	format. 

	Input file must be .fasta or .phylip alignment sequences. The function has two
	options to remove stops codons: the first one removes all the sequences 
	with more than one stop codon and the other one removes the position in all the sequences
	where a stop codon is found. By default the function doesn't remove stop codons.
	The function can make a zone selection for this the user introduce
	a 0 and 1 sequence as a plane file text (-binary name.txt). Finally, if the user 
	wants the protein alignment(-protalign), the function will translate the codon alignment
	with the option of choosing a desired codon table or leaving the default one. 
	

	Arguments:
	-seqdir- alignment file 
	-protalign- protein alignment
	-outputname- exit file name
	
	-outputnamealign- protein alignment file name
	-codon_table- genetic code number (http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
	-stops- c or d. "c" removes all the sequences with premature stop codons. "d" removes the columns of the aligment
	that follows stop codons, including that codon."""

def parcial(seqdir,protalign,outputname,tupla,outputnamealign,codon_table,binary,stops):

    """This funtion take the arguments from the bash and send it to the diferent modules
	 

	>>> parcial("./examples/example.fasta",False,"alineamiento.fasta",None,"alineamientoprot.fasta",1,None,None)
	>>> parcial("./examples/example.fasta",True,"alineamiento.fasta",None,"alineamientoprot.fasta",1,None,None)

	>>> parcial("./examples/example.fasta",False,"alineamiento.fasta",None,"alineamientoprot.fasta",1,"./examples/binary.txt",None)
	>>> parcial("./examples/example.fasta",False,"alineamiento.fasta",None,"alineamientoprot.fasta",1,None,"c")
	>>> parcial("./examples/example.fasta",False,"alineamiento.fasta",None,"alineamientoprot.fasta",1,None,"d")"""


#Here we import the modules
    from modules import Infile
    from modules import Selectingbyzones
    from modules import gap_cleaner
    from modules import removestops
    from modules import alignproteins
    from modules import outfile
    from modules import ID
    
    if seqdir == None:
        raise Exception("wrong argument -seqdir: it should be input.fasta or input.phylip, instead it was " + str(seqdir))

    if binary == None or "." in binary:
	pass
    else:
        raise Exception("wrong argument -binary: it should be name.txt , instead it was "+str(binary)) 

    if type(codon_table) == int:
        pass
    else:
        raise Exception("wrong argument -codon_table: it should be 1,2,3,4,5,6,7,8,9,10 or 11, instead it was "+str(codon_table)) 


    if stops == "d" or stops == "c" or stops == None:
        pass
    else:
        raise Exception("wrong argument -stops: it should be c or d, instead it was "+str(stops))   
    ID = ID._ID(seqdir)
    array =Infile._input(seqdir)
#We send it  to the cutting module
    
    nostops = removestops._remove_stops(array,codon_table,ID,stops)
    nostopsarray = nostops[0]
    ID = nostops[1]
#We send it to the gap cleaner module
    interestarray = Selectingbyzones._zoneselector(nostopsarray,tupla,binary)
    
#We send it to the sotp codones cleaner module
    final = gap_cleaner._gap_cleaner(interestarray)
    
#Ir returns the alignment protein, when necessary
    if protalign == True:
        arrayprotalign = alignproteins._Alignproteins(final,codon_table)   

#The user transforms it into the desired format and this is sent to a file
    
    outfile._outfile(final,outputname,ID)

    if protalign == True:
    	outfile._outfile(arrayprotalign,outputnamealign,ID)
               
            
        
        
    
	

if __name__ == "__main__":
    
	 

	
    import argparse
#These are the arguments that the program will take in

    parser = argparse.ArgumentParser(description="This function cleans codon alignments by gaps and stop codons removal, allowing the user to select zones, generate a protein alignment and select output name and format. 	Input file must be .fasta or .phylip alignment sequences. The function has two options to remove stops codons: the first one removes all the sequences	with more than one stop codon and the other one removes the position in all the sequences where a stop codon is found. By default the function doesn't remove stop codons. The function can make a zone selection for this the user introduce a 0 and 1 sequence as a plane file text (-binary name.txt). Finally, if the user wants the protein alignment(-protalign), the function will translate the codon alignment with the option of choosing a desired codon table or leaving the default one.")
	
    parser.add_argument('--seqdir',  type=str, help='input name example: -seqdir input.fasta')
    parser.add_argument('--stops',  type=str, help='type of stops codons treatment "c" removes all the sequences with premature stop codons. "d" removes the columns of the aligment that follows stop codons, including that codon. example: -stops c')
    parser.add_argument('--binarydir',  type=str, help='.txt that contains a 0 and 1 sequence to select zones of interest example: -binary file.txt')
    parser.add_argument('--protalign', action='store_const', const = True, default = False, help='if True will return cleaned protein aligment ')
    parser.add_argument('--outputname', type=str, default = "parcial.fasta", help="output name. example: -outputname output.fasta")
    parser.add_argument('--outputnamealign', type=str, default = "alignment.fasta", help = "output name for the alignment. example: -outputnamealign alignoutput.fasta")
    parser.add_argument('--tupla', type=str, help='list of lists containing the zones of interest, example: [[3,9],[30,333]]')
    parser.add_argument('--codon_table', type=int, default = 1,help="codon table number, for default its 1. example: -codon_table 11")
#This fits the arguments into a variable    
    args = parser.parse_args()
    args.tupla = None
#This sends the arguments to the function "parcial" as usable arguments for the function
    parcial(args.seqdir,args.protalign,args.outputname,args.tupla,args.outputnamealign,args.codon_table,args.binarydir,args.stops)
