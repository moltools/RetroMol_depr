from retromol.antismash import retrieve_antismash_result, parse_antismash_result


job_id = "bacteria-aa37c9a0-98f8-4273-a59b-b71b8fcee671"
ncbi_acc = "AM420293.1"

src = retrieve_antismash_result(job_id, ncbi_acc)
seqs = parse_antismash_result(src)