import sys, os, json
from sistr.sistr_cmd import main



def test_serotyping(fasta_path):
    sys.argv[1:] = ["-i",fasta_path, "test", "--more-results", "--run-mash", "--qc", "-l", "-vvv", "-o", "sistr_test_results"]
    main()
    assert os.path.exists("sistr_test_results.json"), f"Results SISTR file for {fasta_path} not found"
    with open("sistr_test_results.json") as fp:
        sistr_json_results = json.load(fp)[0]
    assert  sistr_json_results['antigenic_formula'] == "58:l,z13,z28:z6"
    assert  sistr_json_results['cgmlst_matching_alleles'] == 330
    assert  sistr_json_results['cgmlst_found_loci'] == 330
    assert  sistr_json_results['serovar_in_list'] == "Y"
    assert  sistr_json_results['serovar'] == "II 58:l,z13,z28:z6"

def test_noserovarlist_file(fasta_path):
    sys.argv[1:] = ["-i",fasta_path, "test", "--more-results", "--run-mash", "--qc", "-l", "no_file_exists.txt", "-vvv", "-o", "sistr_test_results"]
    main() 
    with open("sistr_test_results.json") as fp:
        sistr_json_results = json.load(fp)[0]
    assert 'serovar_in_list' not in sistr_json_results.keys()       