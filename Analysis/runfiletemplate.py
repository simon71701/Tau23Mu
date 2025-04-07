import subprocess

def get_path(i, era, v):

	file = f'/store/user/bsimon/ParkingDoubleMuonLowMass{i}/SkimPhiPi_2024era{era}_v{v}_stream{i}_Mini_v1'



	output = subprocess.check_output(['gfal-ls', 'davs://eos.cms.rcac.purdue.edu:9000'+file], universal_newline=True)
	_, output = output.split('\n')

	output = output.replace('\n', '')

	new_file = file + '/' + output
	return new_file

