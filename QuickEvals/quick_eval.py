import os
import re
import glob
import time

n_parts = 70

ntupes_info_dir = '/mnt/Lustre/agrp/jonathas/btagntuple_analysis_framework/ntuples_info/'
testfile_with_filenames= ntupes_info_dir+'V47_test_data_partial.txt'
short_name='V47v1454'

##### fixed parameters #####
logdir='/mnt/Lustre/agrp/jonathas/logfiles/'
taggerfile = '\'taggers_v1454.txt\''

# get list of files to evaluate

f=open(testfile_with_filenames)
files_to_eval=f.readlines()
f.close()

#divide list into N parts
def partition(lst, n):
    q, r = divmod(len(lst), n)
    indices = [q*i + min(i, r) for i in xrange(n+1)]
    return [lst[indices[i]:indices[i+1]] for i in xrange(n)]

divided_files = partition(files_to_eval,n_parts)

#create N text files listing the files, save them inside QuickEvals/FileLists
if not(os.path.exists('FileLists')):
  print 'creating directory : FileLists'
  os.mkdir('FileLists')
if not(os.path.exists('Parts')):
  print 'creating directory : Parts'
  os.mkdir('Parts')
os.chdir('FileLists')

txt_file_list = []
job_names = []
output_names = []

for i, file_list in enumerate(divided_files):
	txtlist_name = short_name+'_part'+str(i)+'.txt'
	txt_file_list.append(txtlist_name)

	job_names.append(short_name+'_part'+str(i))
	output_names.append(short_name+'_part'+str(i)+'.root')

	f=open(txtlist_name,'w')
	for file_name in file_list:
		f.write(file_name)
	f.close()

os.chdir('..')

#### variables #####
# job_name = 'quickevaltest'
# outputname = '\'quicktest.root\''
# txtlist = '\'quickevalTrial.txt\''

for i, txtlist in enumerate(txt_file_list):
	job_name = job_names[i]
	outputname = output_names[i]

	bashcommand='qsub -o '+logdir+job_name+'_o.txt -e '+logdir+job_name+'_e.txt -q S -N '+job_name+' evalscript.sh -v OUTPUTNAME='+outputname+',TAGGERFILE='+taggerfile+',TXTLIST='+txtlist
	os.system(bashcommand)
	time.sleep(3)

