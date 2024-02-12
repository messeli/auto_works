# This file helps with finding the labels on the 'UZ' samplings in AUTO
# First the bifurcation diagram b.xxx is opened and regex is used 
# to find the lines of mathing Omeg. All of the lines are printed in a 
# file. Then the pattern for specific amplitude research is returned.
#
# For SublimeText3 Python3 building, the below works : 
# {
#     "cmd": ["python3", "$file"],
#     "selector": "source.python",
#     "file_regex": "^\\s*File \"(...*?)\", line ([0-9]*)"
# }


import re 
#| Method 1 : with open() as f : ... #as far as I remember this did not work in .auto file. In any case you will run the below .py code separately
#| Method 2 : f = open() , f.close() 
def writetofilematchingOmegpattern_generateAmppattern(Omeg,Amp):
	f = open('b.all_sideStudy2','r')
	content = f.read()
	# f.seek(0)
	# lines = f.readlines()

	pattern = r"\n.{13}4.{7}("+Omeg+").+?\n" #|4:UZ
	
	pat = re.compile(pattern)
	matches = pat.finditer(content)

	# spans = []
	# matchTexts = []
	text = "" 
	for i in matches: 
	   # spans.append(i.span(0))
	   # matchTexts.append(i.group(0))
	   text = text + i.group(0) #|group(0): all pattern match.
	f.close()

	f = open('myDummyWriteFile.txt', "w")
	f.write(text)
	f.close()	

	patternForAmp = r"\n.{13}4.{7}("+Omeg+").+?("+Amp+").+?\n"
	return patternForAmp

print(writetofilematchingOmegpattern_generateAmppattern("3.499","1.25"))
