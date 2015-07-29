#! /usr/bin/python
import optparse
import sys
import logging
import time
from util.utils import *
from util.info_ext import *
from mapping.bwa import *
from util.sample import *
from mapping.picard import *
from mapping.gatk import *
from mutation.somatic import *
from SV.breakmer import *
from annotation.annotation import *

usage = 'usage: %prog [options] arg1 arg2'
parser = optparse.OptionParser(usage=usage, version='%prog v3.0 made by coonya')
parser.add_option('-d', dest='basecall',help='basecall DIR')
(options, args) = parser.parse_args()



def check_option(option, msg):
	if not (option):
		print msg
		sys.exit(parser.print_help())

def bwa_aln(sample):
	bwa.aln1(sample)
	bwa.aln2(sample)

def mapping(sample):
	cmd = []

	cmd.append(bwa.aln1(sample))
	cmd.append(bwa.aln2(sample))
	
	exec_command = '\n'.join(cmd)
	
	return exec_command

def postmapping(sample):

	cmd = []

	cmd.append(bwa.sampe(x))
	cmd.append(picard.MergeBamAlignment(x))
	cmd.append(picard.dedup(x))
	cmd.append(gatk.IndelRealigner(x))
	cmd.append(gatk.CountCovariates(x))
	cmd.append(gatk.TableRecalibration(x))

	cmd2 = '\n'.join(cmd)

	return cmd2

def aggregation(tumor, normal=''):
	cmd = []

	if normal != '':
		cmd.append(gatk.RealignerTargetCreator(tumor, normal))
		cmd.append(gatk.IndelRealigner_agg(tumor, normal))
		#cmd.append(picard.CalculateHsMetrics(tumor, 'exon'))
		#cmd.append(picard.CalculateHsMetrics(normal, 'exon'))
		#cmd.append(picard.CalculateHsMetrics(tumor, 'full'))
		#cmd.append(picard.CalculateHsMetrics(normal, 'full'))
		#cmd.append(picard.CalculateHsMetrics(tumor, 'intron'))
		#cmd.append(picard.CalculateHsMetrics(normal, 'intron'))
	else:
		cmd.append(copy_unmateched_normal(tumor, working_dir))
		cmd.append(gatk.RealignerTargetCreator(tumor))
		cmd.append(gatk.IndelRealigner_agg(tumor))
		#cmd.append(picard.CalculateHsMetrics(tumor, 'exon'))
		#cmd.append(picard.CalculateHsMetrics(tumor, 'full'))
		#cmd.append(picard.CalculateHsMetrics(tumor, 'intron'))
	
	exec_command = '\n'.join(cmd)

	return exec_command



def calculate_metrics(sample):
	cmd = []

	cmd.append(picard.CalculateHsMetrics(sample, 'exon'))
	cmd.append(picard.CalculateHsMetrics(sample, 'full'))
	cmd.append(picard.CalculateHsMetrics(sample, 'intron'))

	exec_command = '\n'.join(cmd)

	return exec_command

	

def callingmutation(tumor, normal=''):
	cmd = []

	if normal != '':
		cmd.append(somatic.mutect(tumor, normal))
		cmd.append(somatic.SI(tumor, normal))
		cmd.append(anno.vcf2maf(tumor))
		cmd.append(anno.maf_filter(tumor))

	else:
		cmd.append(somatic.mutect(tumor))
		cmd.append(somatic.SI(tumor))
		cmd.append(anno.vcf2maf(tumor))
		cmd.append(anno.maf_filter(tumor))

	exec_command = '\n'.join(cmd)

	return exec_command

def get_current_job_id(cmd):
	splited_cmd = cmd.split(' ')
	p = subprocess.Popen(splited_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout,stderr = p.communicate()
	stdout = stdout.split('\n')[0].split('.')[0]

	return stdout

	
cpu_no = {
	'basecall': 16,
	'SamToFastq': 2,
	'aln': 16,
	'postmapping': 2,
	'agg': 1,
	'callingmutation': 1,
	'SV': 8,
	'calculate_metrics': 1,
	'report': 1}



########## path or file ################
## config file
config_file = '/home/D131748/src/pipeline/pipeline_config.cfg'

## path for pipeline
working_dir = '/home/D131748/Research/OncoPanel/lane_level/%s' % options.basecall
log_dir = '%s/log' % working_dir
config_dir = '%s/configure' % working_dir
config_dir = '%s/configure' % working_dir
mutect_dir = '%s/mutect' % working_dir
SI_dir = '%s/SI' % working_dir
metrics_dir = '%s/metrics' % working_dir
sample_file = '%s/sample.txt' % config_dir
agg_dir = '%s/tmp/aggregation' % working_dir

logo()



########## checking options or directories or files #####################
check_option(options.basecall, "You need to use '-d' option with basecall DIR\n")
check_working_dir(working_dir)
check_dir(log_dir)
check_dir(mutect_dir)
check_dir(SI_dir)
check_dir(metrics_dir)
check_dir(agg_dir)



################# Generation of logger ################
os.chdir(log_dir)
setup_logger(log_dir, 'root')
logger = logging.getLogger('root')



############# extract parameters ###################
params = read_config(config_file)
logger.info('A extraction of parameters is done.')



############### extract information ##################
info_ext = extract(options.basecall).ext_info()
logger.info('A extraction of information is done.')



################# R1 : get sample information ########
(sample_dic, sample_list, normal_list, tumor_list, paired_list, unpaired_list) = get_sample_list(sample_file)
print sample_dic
file_dic = get_file_name(options.basecall, sample_list)
print_sample(sample_list, normal_list, tumor_list, paired_list, unpaired_list)



############## Generation of picard object
picard = picard(params, options.basecall, file_dic)
bwa = bwa(params, file_dic)
gatk = gatk(params, options.basecall, file_dic)
somatic = somatic(params, file_dic)
breakmer = breakmer(options.basecall, file_dic)
anno = annotation(params, file_dic)



################ R2 : Extract Illumina Barcodes and Illumina Basecalls To Sam ####################
job_id = get_start_jobid()

cmd = []
cmd.append(picard.ExtractIlluminaBarcodes())
cmd.append(picard.IlluminaBasecallsToSam())
cmd2 = '\n'.join(cmd)

file_name = '01.BaseCall.%s' % job_id
wfile = open('%s/%s.sh' % (log_dir, file_name),'w')
wfile.write('%s\n' % pbs_header1(file_name, cpu_no['basecall'], cmd2))
wfile.close()

cmd3 = 'qsub %s/%s.sh' % (log_dir, file_name)
job_id = get_current_job_id(cmd3)

qsub_key = 'basecall'
qsub_id = {qsub_key : job_id}

print_info('Submitting job for Illumina Basecalls to SAM is done.')


sample_list.sort()


for x in sample_list:

	########### Sam to Fastq #####################
	job_id = get_start_jobid()
	wait_id = qsub_id['basecall']

	file_name = '02.SamToFastq.%s.%s' % (x, job_id)
	wfile = open('%s/%s.sh' % (log_dir, file_name),'w')
	wfile.write('%s\n' % pbs_header2(file_name, cpu_no['SamToFastq'], 'no', wait_id, picard.SamToFastq(x)))
	wfile.close()

	cmd3 = 'qsub %s/%s.sh' % (log_dir, file_name)
	job_id = get_current_job_id(cmd3)

	qsub_key = 'SamToFastq_%s' % x
	qsub_id[qsub_key] = job_id


	########### bwa mapping ########################
	job_id = get_start_jobid()
	wait_id = qsub_id[qsub_key]

	file_name = '03.bwa_aln.%s.%s' % (x, job_id)
	wfile = open('%s/%s.sh' % (log_dir, file_name),'w')
	wfile.write('%s\n' % pbs_header2(file_name, cpu_no['aln'], 'no', wait_id, mapping(x)))
	wfile.close()

	cmd3 = 'qsub %s/%s.sh' % (log_dir, file_name)
	job_id = get_current_job_id(cmd3)

	qsub_key = 'bwa_aln_%s' % x
	qsub_id[qsub_key] = job_id



	############ post mapping procedure #############
	job_id = get_start_jobid()
	wait_id = qsub_id[qsub_key]

	file_name = '04.post_mapping.%s.%s' % (x, job_id)
	wfile = open('%s/%s.sh' % (log_dir, file_name),'w')
	wfile.write('%s\n' % pbs_header2(file_name, cpu_no['postmapping'], 'no', wait_id, postmapping(x)))
	wfile.close()

	cmd3 = 'qsub %s/%s.sh' % (log_dir, file_name)
	job_id = get_current_job_id(cmd3)

	qsub_key = 'postmapping_%s' % x
	qsub_id[qsub_key] = job_id

	
	print_info('Submitted job for SamToFastq, BWA mapping and post mapping of %s.' % x)



############# aggregation ###########################
jid_agg = []

for x in sample_dic:
	if len(sample_dic[x]) == 2:
		tumor = sample_dic[x]['TUMOR']
		normal = sample_dic[x]['NORMAL']
		sample1 = sample_dic[x]['TUMOR'].split('_')
		sample = '%s_%s' % (sample1[0], sample1[2][:-1])

		job_id = get_start_jobid()
		qsub_key1 = 'postmapping_%s' % tumor
		qsub_key2 = 'postmapping_%s' % normal
		wait_id = '%s,%s' % (qsub_id[qsub_key1], qsub_id[qsub_key2])
	
		file_name = '05.aggregation.%s.%s' % (sample, job_id)
		wfile = open('%s/%s.sh' % (log_dir, file_name),'w')
		wfile.write('%s\n' % pbs_header2(file_name, cpu_no['agg'], 'yes', wait_id, aggregation(tumor, normal)))
		wfile.close()
	
		cmd3 = 'qsub %s/%s.sh' % (log_dir, file_name)
		job_id = get_current_job_id(cmd3)
	
		qsub_key = 'agg_%s' % tumor
		qsub_id[qsub_key] = job_id
		jid_agg.append(str(job_id))

		print_info('Submitted job for aggregation of %s.' % sample)

	elif len(sample_dic[x]) == 1 and sample_dic[x].keys()[0] == 'TUMOR':
		tumor = sample_dic[x]['TUMOR']
		sample1 = sample_dic[x]['TUMOR'].split('_')
		sample = '%s_%s' % (sample1[0], sample1[2][:-1])

		job_id = get_start_jobid()
		qsub_key1 = 'postmapping_%s' % tumor
		wait_id = qsub_id[qsub_key1]
	
		file_name = '05.aggregation.%s.%s' % (sample, job_id)
		wfile = open('%s/%s.sh' % (log_dir, file_name),'w')
		wfile.write('%s\n' % pbs_header2(file_name, cpu_no['agg'], 'no', wait_id, aggregation(tumor)))
		wfile.close()
	
		cmd3 = 'qsub %s/%s.sh' % (log_dir, file_name)
		job_id = get_current_job_id(cmd3)

		qsub_key = 'agg_%s' % tumor
		qsub_id[qsub_key] = job_id
		jid_agg.append(str(job_id))

		print_info('Submitted job for aggregation of %s.' % sample)


################ calling mutation #####################
for x in sample_dic:
	if len(sample_dic[x]) == 2:
		tumor = sample_dic[x]['TUMOR']
		normal = sample_dic[x]['NORMAL']
		sample1 = sample_dic[x]['TUMOR'].split('_')
		sample = '%s_%s' % (sample1[0], sample1[2][:-1])

		job_id = get_start_jobid()
		qsub_key1 = 'agg_%s' % tumor
		wait_id = qsub_id[qsub_key1]
	
		file_name = '06.CallingMutation.%s.%s' % (sample, job_id)
		wfile = open('%s/%s.sh' % (log_dir, file_name),'w')
		wfile.write('%s\n' % pbs_header2(file_name, cpu_no['callingmutation'], 'no', wait_id, callingmutation(tumor, normal)))
		wfile.close()
	
		cmd3 = 'qsub %s/%s.sh' % (log_dir, file)
		job_id = get_current_job_id(cmd3)
	
		qsub_key = 'callingmutation_%s' % tumor
		qsub_id[qsub_key] = job_id

		print_info('Submitted job for calling mutations of %s.' % sample)

	elif len(sample_dic[x]) == 1 and sample_dic[x].keys()[0] == 'TUMOR':
		tumor = sample_dic[x]['TUMOR']
		sample1 = sample_dic[x]['TUMOR'].split('_')
		sample = '%s_%s' % (sample1[0], sample1[2][:-1])

		job_id = get_start_jobid()
		qsub_key1 = 'agg_%s' % tumor
		wait_id = qsub_id[qsub_key1]
	
		file_name = '06.CallingMutation.%s.%s' % (sample, job_id)
		wfile = open('%s/%s.sh' % (log_dir, file_name),'w')
		wfile.write('%s\n' % pbs_header2(file_name, cpu_no['callingmutation'], 'no', wait_id, callingmutation(tumor)))
		wfile.close()
	
		cmd3 = 'qsub %s/%s.sh' % (log_dir, file_name)
		job_id = get_current_job_id(cmd3)

		qsub_key = 'callingmutation_%s' % tumor
		qsub_id[qsub_key] = job_id

		print_info('Submitted job for calling mutations of %s.' % sample)


###### calling SV #################
qsub_id_all = ','.join(jid_agg)
job_id = get_start_jobid()
file_name = '07.01.callingSV_prep.%s' % job_id
wfile = open('%s/%s.sh' % (log_dir, file_name),'w')
wfile.write('%s\n' % pbs_header2(file_name, cpu_no['SV'], 'yes', qsub_id_all, breakmer.prep()))
wfile.close()

cmd3 = 'qsub %s/%s.sh' % (log_dir, file_name)
job_id = get_current_job_id(cmd3)

qsub_key = 'callingSV_prep'
qsub_id[qsub_key] = job_id


for z in sample_list:
	job_id = get_start_jobid()
	file_name = '07.02.callingSV.%s.%s' % (z, job_id)
	wfile = open('%s/%s.sh' % (log_dir, file_name),'w')
	wfile.write('%s\n' % pbs_header2(file_name, cpu_no['SV'], 'no', qsub_id['callingSV_prep'], breakmer.run_breakmer(z)))
	wfile.close()
	
	cmd3 = 'qsub %s/%s.sh' % (log_dir, file_name)
	job_id = get_current_job_id(cmd3)

	qsub_key = 'SV_%s' % z
	qsub_id[qsub_key] = job_id

	print_info('Submitted job for calling SV of %s.' % z)

###### calculate metrics #################

jid_cal_metrics = []

for z in sample_list:
	job_id = get_start_jobid()

	file_name = '08.calculate_metrics.%s.%s' % (log_dir, file_name)
	wfile = open('%s/%s.sh' % (log_dir, file_name),'w')
	qsub_key = 'SV_%s' % z
	wfile.write('%s\n' % pbs_header2(file_name, cpu_no['calculate_metrics'], 'no', qsub_id[qsub_key], calculate_metrics(z)))
	wfile.close()
	
	cmd3 = 'qsub %s/%s.sh' % (log_dir, file_name)
	job_id = get_current_job_id(cmd3)

	qsub_key = 'calculate_metrics_%s' % z
	qsub_id[qsub_key] = job_id
	jid_cal_metrics.append(str(job_id))

	print_info('Submitted job for calculating metrics of %s.' % z)

###### Generating excel report #################
qsub_id_all = ','.join(jid_cal_metrics)

file_name = '09.report.%s' % options.basecall.replace('/','')
wfile = open('%s/%s.sh' % (log_dir, file_name),'w')
wfile.write('%s\n' % pbs_header2(file_name, cpu_no['report'], 'yes', qsub_id_all, excel_report(working_dir)))
wfile.close()

cmd3 = 'qsub %s/%s.sh' % (log_dir, file_name)
job_id = get_current_job_id(cmd3)

print_info('Submitted job for Generating excel report.')
print_info('Complete to submit jobs to AMC clusters.')
