#!/usr/bin/env python3

import click
import json
import os
import re
import shutil

# Ensure files adhere to the rules defined in the the chewbbaca allele call function

def chewbbaca_filename_format(filename):
    # Rule 1 Replace spaces and illegal characters with underscores
    name, ext = os.path.splitext(filename)
    new_filename = re.sub(r'[^A-Za-z0-9_\-]', '_', name)
    new_filename = "{}{}".format(new_filename, ext)
    return new_filename

def copy_new_file(clean_fasta_dir, new_name, filename, original_path):
    # deal with moving the files 
    clean_path = os.path.join(clean_fasta_dir, new_name)
    # If the filename was changed, copy the renamed file to the output directory
    if filename != new_name:
        print("Renaming and copying: {} -> {}".format(filename, new_name))
        shutil.copy2(original_path, clean_path)
    else:
        print("Copying: {}".format(filename))
        shutil.copy2(original_path, clean_path)

@click.group()
def cli():
    """ This script supports the Core Genome MLST service with multiple commands."""
    pass

@cli.command()
@click.argument("service_config")
def clean_fasta_filenames(service_config):
    """Ensure files adhere to the rules defined by chewbbacca"""

    with open(service_config) as file:
        data = json.load(file)
        raw_fasta_dir = data["raw_fasta_dir"]
        clean_fasta_dir = data["clean_data_dir"]
        for filename in sorted(os.listdir(raw_fasta_dir)):
            original_path = os.path.join(raw_fasta_dir, filename)
            new_name = chewbbaca_filename_format(filename)
            copy_new_file(clean_fasta_dir, new_name, filename, original_path)

if __name__ == "__main__":
    cli()