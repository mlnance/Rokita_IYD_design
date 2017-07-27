import numpy as np
import matplotlib.pyplot as plt
import sys, os

def get_decoy_list(dir):

	decoy_list = [os.path.join(dir, "starting_conf.pdb.ppk")]

	for decoy in sorted( os.listdir(dir) ):
		if decoy.endswith(".ppk") and decoy != "starting_conf.pdb.ppk" and not decoy.startswith("._"):
			decoy_list.append( os.path.join(dir, decoy) )

	return decoy_list

def get_rmsd(decoy, bound, partner):

	if partner == 1:
		intf = " and chain A and resi 8+32+36+39+40+42+45+49+106+107+108+137+138+139+140+143+190+192+193+194+198+221+222+223+224+225+228+246+248+249+250+274 and name CA+C+N"

	elif partner == 2:
		intf = " and chain D and resi 36+37+40+41+43+46+50+110+111+112+113+147+148+151+152+196+197+198+199+202+229+231+232+233+234+238+239+256+258+259+260+261+262+284 and name CA+C+N"

	else:
		sys.exit("Unreachable code. Problem with partner number.")
		
	rmsd = cmd.align(decoy + " and name CA+C+N", bound + " and name CA+C+N", cycles = 0)[0]
	intf_rmsd = cmd.align(decoy + intf, bound + intf, cycles = 0)[0]

	return rmsd, intf_rmsd

def rmsd_to_bound(decoy_list, bound, partner):

	decoy_stats = []
	decoy_number = 0
	#decoy_stats_file = open("decoy_stats", 'w')

	for decoy in decoy_list:
		decoy_number += 1
		cmd.load(decoy)
		decoy_name = decoy.split("/")[-1]
		rmsd, intf_rmsd = get_rmsd(decoy_name, bound, partner)
		cmd.delete(decoy_name)
		decoy_stats.append( [decoy_name, decoy_number, rmsd, intf_rmsd] )
		#decoy_stats_file.write(decoy_name + "\t" +  str(decoy_number)  + "\t" + str(rmsd) + "\t" + str(intf_rmsd) + "\n")
		#decoy_stats_file.flush()

	#decoy_stats_file.close()
	return decoy_stats

def map_num_to_rmsd(num, decoy_stats):

	for decoy in decoy_stats:
		if num == decoy[1]:
			return decoy[2], decoy[3]

	sys.exit("Decoy" + str(num) + " not found!")


def parse_score_file(filename, n, rec_stats, lig_stats):

	print "Parsing file: " + filename

	scorefile = open(filename, 'r').readlines()
	incorrect = []
	acceptable = []
	medium = []
	high = []

	if n == -1:
		n = len(scorefile) - 2
	elif n > len(scorefile) - 2:
		print "The score file only has " + str (len(scorefile) - 2) + "entries."
		n = len(scorefile) - 2 

	for i in range( 2, n+2 ):
		score_split = scorefile[i].split()
		
		if len(score_split) == 40 or len(score_split) == 38:
			# makes it easier to determine y_axis limits later

			rec_num = int(float(score_split[10]))
			lig_num = int(float(score_split[11]))

			rec_rmsd, rec_intf_rmsd = map_num_to_rmsd(rec_num, rec_stats)
			lig_rmsd, lig_intf_rmsd = map_num_to_rmsd(lig_num, lig_stats)

			if score_split[3] == "0.000":
				incorrect.append( [ score_split[-1], rec_rmsd, rec_intf_rmsd, lig_rmsd, lig_intf_rmsd ] )
			elif score_split[3] == "1.000":
				acceptable.append( [ score_split[-1], rec_rmsd, rec_intf_rmsd, lig_rmsd, lig_intf_rmsd ] )
			elif score_split[3] == "2.000":
				medium.append( [ score_split[-1], rec_rmsd, rec_intf_rmsd, lig_rmsd, lig_intf_rmsd ] )
			elif score_split[3] == "3.000":
				high.append( [ score_split[-1], rec_rmsd, rec_intf_rmsd, lig_rmsd, lig_intf_rmsd ] )
			else:
				sys.exit("unreachable capri rank")

		else:
			print "There is an issue with line number: " + str(i+1)

	incorrect = map( list, zip(*incorrect) )
	if len(acceptable) > 0:
		acceptable = map( list, zip(*acceptable) )
	if len(medium) > 0:
		medium = map( list, zip(*medium) )
	if len(high) > 0:
		high = map( list, zip(*high) )

	scores = [incorrect, acceptable, medium, high]

	return scores

def make_plots(scores):

	plt.figure(1, figsize = (8,8))

	if len(scores[0]) > 0:
		incorrect = plt.scatter( scores[0][2], scores[0][4], color=(0.75, 0.75, 0.75), alpha=0.8, edgecolors=(0.2, 0.2, 0.2) )
	if len(scores[1]) > 0:
		acceptable = plt.scatter( scores[1][2], scores[1][4], color=(1.00, 0.75, 0.25), alpha=0.8, edgecolors=(0.2, 0.2, 0.2) )
	if len(scores[2]) > 0:
		medium = plt.scatter( scores[2][2], scores[2][4], color=(1.00, 0.25, 0.00), alpha=0.8, edgecolors=(0.2, 0.2, 0.2) )
	if len(scores[3]) > 0:
		high = plt.scatter( scores[3][2], scores[3][4], color=(0.00, 0.75, 0.00), alpha=0.8, edgecolors=(0.2, 0.2, 0.2) )

	plt.title("Target: 2J7P")

	plt.xlim( [1.0, 5.0] )
	plt.ylim( [1.0, 5.5] )

	plt.xlabel("Receptor RMSD" + r'$_c$' + " (in " + r'$\AA$' + ")", fontsize=16)
	plt.ylabel("Ligand RMSD" + r'$_c$' + " (in " + r'$\AA$' + ")", fontsize=16)

	var_list = []
	label_list = []

	if 'incorrect' in locals():
		var_list.append(incorrect)
		label_list.append("Incorrect")

	if 'acceptable' in locals():
		var_list.append(acceptable)
		label_list.append("Acceptable")

	if 'medium' in locals():
		var_list.append(medium)
		label_list.append("Medium quality")

	if 'high' in locals():
		var_list.append(high)
		label_list.append("High quality")

	plt.legend( var_list, label_list, scatterpoints=1, loc=1, ncol=4, fontsize=15, fancybox=True, framealpha=0.5 )
	plt.savefig( "2J7P_doped_RMSD_to_bound", bbox_inches="tight", transparent=True, dpi=300 )


cmd.load("2J7P_b.pdb")
bound = "2J7P_b"
rec_dir = "rec_ensemble"
lig_dir = "lig_ensemble"

rec_list = get_decoy_list(rec_dir)
lig_list = get_decoy_list(lig_dir)

rec_stats = rmsd_to_bound(rec_list, bound, 1)
lig_stats = rmsd_to_bound(lig_list, bound , 2)

scores = parse_score_file("docking_score.sc", 5000, rec_stats, lig_stats)
make_plots(scores)
