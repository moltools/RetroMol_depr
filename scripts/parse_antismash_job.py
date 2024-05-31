import argparse 
from retromol.npkg.antismash import parse_antismash_json, get_antismash_data

parser = argparse.ArgumentParser()
parser.add_argument("--job", required=True, type=str)
args = parser.parse_args()
job_id = args.job
data = get_antismash_data(job_id)
result = parse_antismash_json(data)
print(result)
