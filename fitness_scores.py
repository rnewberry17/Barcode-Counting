# By Robert W. Newberry, PhD
# Last Update: 190909 15:46
#
# Purpose: calculate and plot fitness scores from deep mutational scanning time course experiments
#
# Inputs
#    1) A pickle dictionary relating barcodes to protein variants
#	Keys should be 20bp barcodes, though the length can be changed in the count_barcodes function
#	Values should be tuples, the first element being the position in sequence,
#	    the second being the amino acid in one letter code ('X' for stop)
#	    **WT is represented as (0,'WT')
#    3) Lists of barcodes obtained by next-gen sequencing, each on a separate line
#	    **Quality filtering should be performed separately to generate the lists
#	    **This script expects evenly spaced time points
#		which can be changed by adjusting the calc_fitness function
#	    **Files must be entered in order (i.e. first time point first, etc.)
#	    **This script expects the barcode sequence to be the first 20bp of each read
#		which can be changed by adjusting the count_barcodes function
#    4) A pickle dictionary relating amino acids to positions in a list,
#	    which is used to order the amino acids on the resulting heatmap
#	    **Must be named 'aminotonumber.pkl'
#    5) A pickle dictionary relating amino acids to positions in the protein,
#	    which is used to color the WT residues in the heatmap
#	    **Must be named 'wt_residues.pkl'
#
# Outputs
#    1) A text file containing the barcode, associated variant, and count at each time point,
#	    titled 'merged_barcode_counts.txt'
#    2) A text file containing the fitness scores for each mutant, averaged over all barcodes
#	    mapping to that mutant, titled 'fitness_scores.txt'
#    3) A pickle file of the fitness scores from output #2 (fitness_scores.pkl)
#    4) A heatmap describing the fitness scores for each mutation at each position
#
# Usage
#    1) For each replicate:
#	    a) Run count_barcodes('consensus_dictionary.pkl',*barcode_lists) # substitute filenames as appropriate
#	    b) rename output files
#    2) Create a text file, on each line, include the following
#	    merged_barcode_counts.txt,read_counts.txt # substitute filenames as apporpriate
#    3) Run calc_fitness(Filename_from_step_2,count_threshold) # substitute filename as appropriate
#	    The count_threshold excludes barcodes with insufficient representation in the naive library
#	    We typically use a threshold of 10
#    4) Run make_heatmap(AA_number,WT_median)
#	    where AA_number is the number of residues
#	    and WT_median is returned by calc_fitness in step 3


import cPickle as pickle
import math
from scipy import stats
import statistics
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid


# Count the number of times each barcode appears in each file
def count_barcodes(dictionary,*filenames): # read in a barcode dictionary and 3 barcode lists from counting experiments
    print 'Counting barcodes; this may take some time'
    barcode_dictionary=pickle.load(open(dictionary,'r'))
    with open('merged_barcode_counts.txt','w') as fout, open('read_counts.txt','w') as fout2:

	read_counts=[] # create an empty list for the total number of reads at each time point
        filenumber=0 # for iterating over each of the lists of barcodes
        barcode_counts={} # create an empty dictionary for the barcode count arrays

        for key in barcode_dictionary: # for every barcode in the dictionary
            barcode_counts[key]=[] # create an empty array for the counts to go
	    for x in range(len(filenames)): # for each timepoint
		barcode_counts[key].append(0) # set the initial counts to 0

        for filename in filenames: # for each list of barcodes from a counting experiment
	    counter=0
            with open(filename,'r') as barcode_reads: 
                for line in barcode_reads: # for each barcode
                    counter+=1 # increment the counter
                    try: 
                        barcode_counts[line[0:20]][filenumber]+=1 # increment the count of that barcode in the count dictionary
                    except KeyError: # unless that barcode isn't in the dictionary
                        donothing=0 # in which case, ignore it
	    fout2.write(str(counter)+'\n') # output the total number of reads for the timepoint
            filenumber+=1 # go to the next list of barcodes

        for key in barcode_counts: # for each barcode in the dictionary
            fout.write(key+','+str(barcode_dictionary[key])) # write the barcode and variant to a text file
            for x in range(len(filenames)): # for each time point
                fout.write(','+str(barcode_counts[key][x])) # append the number of times that barcode appears
            fout.write('\n')

    print 'Counting barcodes: complete!'

# Calculate a fitness score for each protein variant
def calc_fitness(filelist_filename, count_threshold): # read in the minimum number of initial counts necssary to include a given barcode in the analysis
    with open(filelist_filename,'r') as filelist:

	barcode_slopes={} # a dictionary to hold the slopes for each barcode in each replicate
	replicate_number=0

	for line in filelist: # for each replicate
	    replicate_number+=1
	    names=line.split(',')
	    barcode_counts_filename=names[0]
	    total_counts_filename=names[1][:-1]

	    with open(barcode_counts_filename,'r') as barcode_counts_file, open(total_counts_filename,'r') as total_counts_file:
	        timepoints=0 # for counting the number of timepoints
	        total_counts=[] # create an empty array for the total counts
	        for line in total_counts_file: # for each line in the total counts file
	            total_counts.append(float(line)) # assign that value to an element of the empty array
	            timepoints+=1 # increase the count of timepoints
	
	        x_values=[] # time points against which to calculate slope
	        for t in range(timepoints):
	            x_values.append(t)

	        count_dict={} # a dictionary to contain the counts for each barcode
	        freq_dict={} # a dictionary to contain the frequency of each barcode at each time point

	        for line in barcode_counts_file: # for each barcode
	            stuff=line.split(',') # put the data for each barcode into an array
	            count_dict[stuff[0]]=[int(stuff[1][1:]),stuff[2][2:-2]] # count_dict[barcode] = [position,AA]
	            for t in range(3,len(stuff)): # then for each timepoint
	                count_dict[stuff[0]].append(float(stuff[t])) # count_dict[barcode] = [position,AA,count1,count2,count3,etc...]
	
	        for barcode in count_dict: # for each barcode
	            if count_dict[barcode][2] > int(count_threshold): # consider only barcodes that start with at least a given number of counts
	                freq_dict[barcode]=[] # create an empty list for the frequencies at each timepoint
	                y_values=[] # create an empty list for the log transformed frequencies that will be regressed against the time course
	                for t in range(2,len(count_dict[barcode])): # for each timepoint
	                    if count_dict[barcode][t] == 0: # if the barcode falls out
	                        count_dict[barcode][t]+=1 # add a pseudocount to allow for log transform
	                    freq_dict[barcode].append(float(count_dict[barcode][t])/total_counts[t-2]) # convert each count to frequency
	                    y_values.append(math.log10(freq_dict[barcode][t-2]))
	                linear_regression=[] # create an empty list to receive the results of the linear regression
	                linear_regression=stats.linregress(x_values,y_values)
			try:
			    barcode_slopes[barcode].append(linear_regression[0])
			except KeyError:
	                    barcode_slopes[barcode]=[count_dict[barcode][0],count_dict[barcode][1],linear_regression[0]] # barcode_slopes[barcode] = [position,AA,slopes...]

	
	barcode_averages={} # a dictionary to contain the average and standard devation of the fitness of a barcode across replicates

	for barcode in barcode_slopes:
	    barcode_averages[barcode]=barcode_slopes[barcode][0:2] 
	    barcode_averages[barcode].append(statistics.mean(barcode_slopes[barcode][2:]))
	    barcode_averages[barcode].append(np.std(barcode_slopes[barcode][2:])) # barcode_averages[barcode] = [position,AA,slope_average,slope_sd]

        set_dict={} # a dictionary to contain all of the slopes for each variant
	fitness_scores={} # a dictionary to contain the fitness scores

        for barcode in barcode_averages: # for each barcode with a slope
            variant=(barcode_averages[barcode][0],barcode_averages[barcode][1]) # variant = [position,AA]
            try: # if there is already a slope for that variant
                set_dict[variant].append((barcode_averages[barcode][2],barcode_averages[barcode][3])) # append the slope for that barcode to the array for slopes for that variant
            except KeyError: # if there is no slope for that variant
                set_dict[variant]=[(barcode_averages[barcode][2],barcode_averages[barcode][3])] # create a new entry in the set dictionary for that variant and append the slope for that barcode

        for variant in set_dict: # for each variant
	    scores=[]
	    errors=[]
	    for entry in set_dict[variant]:
		scores.append(entry[0])
		errors.append((entry[1])**2)
            average=statistics.mean(scores) # calculate the average slope among all of the barcodes mapping to that variant
	    if replicate_number > 1:
		stdev=math.sqrt(sum(errors))/len(errors) # as well as the standard deviation
	    elif replicate_number == 1:
		stdev=np.std(scores)
	    fitness_scores[variant]=[average,stdev,len(scores)] # dump that data into a dictionary
	
# calculate p-values for t-tests compared to WT
	m2=fitness_scores[(0, 'WT')][0]
	s2=fitness_scores[(0, 'WT')][1]
	n2=fitness_scores[(0, 'WT')][2]
	for variant in fitness_scores:
	    m1=fitness_scores[variant][0]
	    s1=fitness_scores[variant][1]
	    n1=fitness_scores[variant][2]
	    try:
		fitness_scores[variant].append(stats.ttest_ind_from_stats(m1,s1,n1,m2,s2,n2,equal_var=False)[1])
	    except ZeroDivisionError:
		fitness_scores[variant].append(np.nan)

	with open('fitness_scores.pkl','w') as fout: # output the data as a pickle
	    fout.write(pickle.dumps(fitness_scores))
	with open('fitness_scores.txt','w') as fout2: # and as text
	    for variant in fitness_scores:
		fout2.write(str(variant)+','+str(fitness_scores[variant])+'\n')

	WT_scores=[]
	for score in set_dict[(0, 'WT')]:
	    WT_scores.append(score[0])
	WT_median=np.median(WT_scores) # calculate the median fitness score of the WT barcodes
	return WT_median


# Make a heatmap for the deep mutational scanning data
def make_heatmap(residue_number,WT_median): # requires the number of residues from user and the WT median from calc_fitness

    fitness_scores=pickle.load(open('fitness_scores.pkl','r')) # load the fitness score file
    aa_to_number=pickle.load(open('aminotonumber.pkl','r')) # load the amino acid order
    wt_residues=pickle.load(open('wt_residues.pkl','r')) # load the wild-type residues

    amino_acids=20
    fitness_array=np.zeros(shape=(amino_acids,residue_number)) # create an array the size of the protein
    fitness_array[:]=np.nan # initalize all the values to NaN

    for variant in fitness_scores: # for each variant
      if variant[1] != 'WT': # for variants other than WT
	fitness_array[aa_to_number[variant[1]]][variant[0]-1]=fitness_scores[variant][0]-WT_median # subtract the WT median
    for position in wt_residues: # for the WT residues
	fitness_array[aa_to_number[wt_residues[position]]][position-1]=0 # set the normalized fitness to 0

    cmap=matplotlib.cm.RdBu # start with a red-white-blue colormap
    color_endpoint=max(abs(np.nanmin(fitness_array)),np.nanmax(fitness_array)) # find the max value to scale the colorbar
    cmap.set_bad(color='lightgrey') # use white for missing values
    masked_array=np.ma.masked_invalid(fitness_array) # mask the missing values
    plt.imshow(masked_array,cmap=cmap, interpolation='nearest') # draw the heatmap
    plt.clim(color_endpoint,color_endpoint) # set the colorbar range to reflect the most extreme value in the data
    plt.colorbar() # include the colorbar
    plt.show() # and show


