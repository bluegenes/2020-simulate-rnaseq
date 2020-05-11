import sys
import argparse

## read in csv, get list of coding, noncoding matches
def summarize_sencha_classification_by_read(sencha_csv, output_summary):
    header = ["read_id", "frames_coding", "frames_stop", "frames_noncoding", "classification", "filename"]
    #initialize vars
    coding_count, noncoding_count, stop_count=0,0,0
    read,classification, origin_file="", "",""

    with open(output_summary, "w") as out:
        out.write("\t".join(header) + "\n")
        # read line
        with open(sencha_csv) as f:
            next(f) #skip header
            for line in f:
                line = line.strip().split(",")
                try:
                    read_id,jaccard_in_peptide_db,n_kmers,category,translation_frame,filename = line #line.strip().split(",")
                except:
                    # some read ids have commas -- pass for now, fix later
                    pass
                    #import pdb;pdb.set_trace()
                # initialize read
                if not read:
                    read=read_id
                    origin_file=filename
               # we've reached a new read
                if read_id != read:
                    if not classification:
                        classification = "non-coding"
                    out.write(read + "\t" + str(coding_count) + "\t" + str(stop_count) + "\t" + str(noncoding_count) + "\t" + classification + "\t" + origin_file + "\n")
                    # reset vars
                    read=read_id
                    origin_file=filename
                    coding_count, noncoding_count, stop_count=0,0,0
                    classification=""
                # now handle this read_id frame
                if category=="Coding":
                   coding_count+=1
                   classification="coding"
                elif category=="Non-coding":
                    noncoding_count+=1
                elif "stop codon" in category:
                    stop_count+=1
                elif "shorter than 3" in category:
                    classification = "NA"
            if not classification:
                classification = "non-coding"
            out.write(read + "\t" + str(coding_count) + "\t" + str(stop_count) + "\t" + str(noncoding_count) + "\t" + classification + "\t" + origin_file + "\n")



if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('csv')
    p.add_argument('-o', '--output')
    args = p.parse_args()
    sys.exit(summarize_sencha_classification_by_read(args.csv, args.output))
