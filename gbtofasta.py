#This python script is very useful when after downloading the sequences from blast for a specific haplogruoup in a genbank file.
#we use this file to retrieve sequence per sequence what is its origin/population and creates a multifasta with this information.
#so the output of running this program corresponds to a multifasta created based on the hits from a sequence selected from our
#dataset belonging to a specific haplogroup of interest, where the first line corresponds to the population/oringin + its original id.
#These resulting files will be added with the multifastas with the sequences from our dataset that belong to that specific 
#haplogroup for the network analysis using Network 10.2.
#as mentioned before, knowing to which population belong each sequence is crucial for the conclusions of this project and
#this script can speed up this process.


import os
from Bio import SeqIO

# Directory path containing GenBank files
directory = "gbdirectory"

# Get the absolute path to the directory
absolute_directory = os.path.abspath(directory)


# Iterate over each GenBank file in the directory
for filename in os.listdir(absolute_directory):
    if filename.endswith(".gb"):
        # Get the full path to the GenBank file
        genbank_file = os.path.join(directory, filename)

        # Create an empty list to store the sequences
        sequences = []

        # Parse the GenBank file
        for record in SeqIO.parse(genbank_file, "genbank"):
            # Retrieve the qualifier (country, isolate, or note) from the features table
            qualifier_value = ""
            qualifiers_to_search = ["country", "isolation_source", "note"]

            for feature in record.features:
                if feature.type == "source":
                    for qualifier in qualifiers_to_search:
                        if qualifier in feature.qualifiers:
                            qualifier_value += feature.qualifiers[qualifier][0] + " "
                            #you get all the qualifiers you find together because you do not know which one corresponds to the population
                            
            # Get the sequence ID
            sequence_id = record.id

            # Get the sequence
            sequence_text = str(record.seq)
			
            if qualifier_value == "":
            # Create the ID for the FASTA entry
                fasta_id = sequence_id
            else:
                fasta_id = qualifier_value + ' ' + sequence_id

            # Create the FASTA entry and add it to the list of sequences
            fasta_entry = f">{fasta_id}\n{sequence_text}"
            sequences.append(fasta_entry)

        # Create the output FASTA filename
        fasta_filename = f"{filename[:-3]}.fasta"  # Remove the extension .gb from the GenBank filename

        # Write the sequences to a FASTA file
        fasta_filepath = os.path.join(directory, "Blast" + fasta_filename)
        with open(fasta_filepath, "w") as fasta_file:
            fasta_file.write("\n".join(sequences))
