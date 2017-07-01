import re
import os
import sys
import glob
import shutil
import subprocess
import datetime
import logging
import logging.handlers


class Blasperin():
    """
    Check configurations, required softwares and run the clustering of EC numbers Fasta files.
    """

    def __init__(self, label=None, author=None, fasta_files_directory=None, log_file=None):

        self.label = label
        self.author = author

        # Keep tracking of what EC files is being clustered.
        self.current_fasta_file = None

        # Fasta source.
        self.fasta_files_directory = fasta_files_directory

        self.log_file = log_file

    def execute_analysis(self):
        """
        Execute the analysis.

        This is the main method that call all other auxiliary clustering methods.
        """

        self.create_log_system(self.log_file)

        self.log.info("-- START --:blasperin:execute_analysis.")

        self.write_metadata()

        self.generate_blast_results()

        self.log.info("-- DONE --:blasperin:execute_analysis.")

    def create_log_system(self, log_file=None):
        """
        Set all logger parameters (file log path, for example), output format and set the class property that stores the logging system.

        """

        log = logging.getLogger('')
        log.setLevel(logging.DEBUG)
        format = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s")

        ch = logging.StreamHandler(sys.stdout)
        ch.setFormatter(format)
        log.addHandler(ch)

        fh = logging.handlers.RotatingFileHandler(
            log_file, maxBytes=0, backupCount=0)
        fh.setFormatter(format)
        log.addHandler(fh)

        self.log = log

    def blast_proteins(self, fasta_file=None):
        """
        Execute BLAST software to generate the similarity results.

        Args:
            fasta_file(str): Path for the Fasta file to be processed.

        """

        if not self.fasta_file_is_done(fasta_file):

            self.purge_temporary_files(fasta_file)

            blast_result_file_name = fasta_file + '.blastall'

            self.log.info("Creating BLAST DB for: " + fasta_file)

            os.popen(
                "makeblastdb -in " +
                fasta_file +
                " -out " +
                fasta_file +
                " -dbtype prot")

            self.log.info("DONE creating BLAST DB for: " + fasta_file)

            self.log.info("BLAST PROTEINS! : " + fasta_file)

            os.popen(
                "blastp -query " +
                fasta_file +
                " -outfmt '6 qseqid sseqid pident ppos score bitscore qstart qend sstart send qlen slen evalue' -evalue 0.1 -num_alignments 10000000 -db " +
                fasta_file +
                " -out " +
                blast_result_file_name)

            self.log.info("DONE BLAST PROTEINS! : " + fasta_file)

            self.write_done_file(fasta_file)

        else:
            self.log.info("File already blasted: " + fasta_file)

    def fasta_file_is_done(self, fasta_file=None):
        """
        Test if the Fasta file was already processed by BLAST.

        Args:
            fasta_file(str): Fasta file path.

        """

        files = self.done_files_list()

        if fasta_file in files:
            return True

    def write_done_file(self, fasta_file=None):
        """
        Update the done_files tracking file.

        Args:
            fasta_file(str): File path to be inserted into the done_files file.
        """

        f = open(self.fasta_files_directory + '/done_files', 'a')

        f.write(fasta_file + "\n")

        f.close()

    def galperin_analysis(self, fasta_file=None):
        """
        Call the BLAST software and call the actual clustering method.

        """

        self.log.info("Start blasting process for: " + fasta_file)

        blast_result = self.blast_proteins(fasta_file)

        self.log.info("End of blasting process for: " + fasta_file)

    def generate_blast_results(self):
        """
        Walk through the directory files and call the 'galpering_analysis' method.

        It also check if the file to be processed was already processed.

        """

        done_files = self.done_files_list()

        files = glob.glob(self.fasta_files_directory + '/' + '*.fasta')

        self.total_of_files = len(files)

        self.log.info("Total of " +
                      str(self.total_of_files) +
                      " will be processed.")

        counter = 1

        for fasta_file in files:

            self.log.info("Processing " + str(counter) + " of " +
                          str(self.total_of_files) + " total files.")
            counter += 1

            fasta_file_name = os.path.basename(fasta_file)

            # Only execute clustering if it wasn't run (and done) before.
            if fasta_file_name not in done_files:
                self.current_fasta_file = os.path.basename(fasta_file)

                self.log.info("Going to cluster: " + fasta_file)

                self.galperin_analysis(fasta_file)
            else:
                self.log.info(
                    "File SKIPED. It was already clustered: " +
                    fasta_file)

    def done_files_list(self):
        """
        Read the 'done_files' and return it files list.

        Returns:
            (list): List of file names.
        """

        finished = []

        done_files = self.fasta_files_directory + '/' + 'done_files'

        if os.path.exists(done_files):
            f = open(done_files)

            for line in f:
                # Remove newline
                line = line.rstrip('\r\n')
                finished.append(line)

            f.close()

        return finished

    def mark_clustering_done(self):
        """
        Write the processed EC number into the 'done_files'.

        This file keep tracking of what was already done.

        """

        done_file = fasta_files_directory + '/' + 'done_files'
        f = open(done_file, 'a')
        f.write(self.current_fasta_file + "\n")
        f.close()

    def purge_temporary_files(self, fasta_file=None):
        """
        Remove temporary BLAST files. Files that ends with .phr, .pin, .psq and .blastall.

        This method is used when a clustering is interrupted for some reason and old processed files

        have to be removed: We have to make sure a whole clustering process was executed neat and clean.
        """

        file_mask = fasta_file

        files = []

        files.append(file_mask + '.phr')
        files.append(file_mask + '.pin')
        files.append(file_mask + '.psq')
        files.append(file_mask + '.blastall')

        for file_to_remove in files:
            if os.path.exists(file_to_remove):
                self.log.info("Purging old file: " + file_to_remove)
                os.remove(file_to_remove)

    def write_metadata(self):
        """
        This is one of the most important methods.

        Metadata file have to be written since that data will be used as a mark for future relational database insertions.

        Without metadata the relational database wouldn't be able to know what clustering method was used.

        """

        metadata_file = self.fasta_files_directory + '/' + 'blasperin_metadata'

        with open(metadata_file, 'w') as f:
            f.write('label = ' + str(self.label) + "\n")
            f.write('author = ' + str(self.author) + "\n")
            f.write('date = ' + str(datetime.datetime.now()) + "\n")
            f.write('software = blast' + "\n")
