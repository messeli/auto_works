import re 
#| Method 1 
# with open('b.all_sideStudy2','r') as f:
# 	content = f.read()
# 	seek(0)
# 	lines = f.readlines()
# 	pattern = r"\n.{12}-4.{7}(7.499).+?\n"
# 	pat = re.compile(pattern)
	
# 	matchs = re.finiter(content)
	
# 	spans = []
#     matchTexts = []
#     text = "" 
#     for i in matches: 
# 	    spans.append(i.span(0))
# 	    matchTexts.append(i.group(0))
# 		text = text + i.group(0)
# with open('myDummyWriteFile.txt', "w") as f:
# 	f.write(text)

#| Method 2
f = open('b.all_sideStudy2','r')
content = f.read()
f.seek(0)
lines = f.readlines()
pattern = r"\n.{13}4.{7}(7.499).+?\n" #|ACTION REQ FOR WHAT TO WRITE
pat = re.compile(pattern)
matches = pat.finditer(content)

spans = []
matchTexts = []
text = "" 
for i in matches: 
   spans.append(i.span(0))
   matchTexts.append(i.group(0))
   text = text + i.group(0) #|Here w/ [:-1], no extra new line
f.close()

f = open('myDummyWriteFile.txt', "w")
f.write(text)
f.close()	