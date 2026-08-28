# Needs to take the S3 location, output save location, and work directory as input arguments

# TODO: Need to have a barcode conversion file!! i.e Ultima to Scale. Probably a .csv file.
# Probably just keep this in the github repo

# Then download the file from S3

# after downloading, then run convertUltimaToIllumina
# will need the path for the downloaded INPUT_FASTQ file, the work directory
# will need to get the PCR barcode from the conversion file.

# after running, then upload the fastq files back to the S3 bucket

# Then clean up the fastq files on the SSD