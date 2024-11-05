while true; do
	case "$1" in
		-i) inputFile=$2; shift 2;;
		-o) outputDir=$2; shift 2;;
		-n) name=$2; shift 2;;
		-a) annopath=$2; shift 2;;
		*) break;
	esac
done


windowSize=1000
stepSize=200
identity=50
queryCover=70
subjectCover=70
threads=40
script=/lustre/home/liudongmei/project/2024/20240625_GCtype3.0/9.circle/cgview/script
database=/data/pipeline/database

	#mkdir -p ${outputDir}/anno
	#prodigal -i ${inputFile} -a ${outputDir}/anno/${name}.prodigal.pep -f gff -o ${outputDir}/anno/${name}.prodigal.gff3
	#/data/public_tools/diamond/bin/diamond blastp -d  /data/pipeline/database/COG/cog.all.dmnd  -q ${outputDir}/anno/${name}.prodigal.pep  -e 1e-5 -p ${threads} --max-target-seqs 1 --id ${identity} --query-cover ${queryCover} --subject-cover ${subjectCover} --sensitive -o ${outputDir}/anno/${name}.out 

	#rnammer -S bac -multi -h ${outputDir}/anno/${name}.RNAmmer.hmmreport -gff ${outputDir}/anno/${name}.RNAmmer.gff -f ${outputDir}/anno/${name}.RNAmmer.fasta ${inputFile}
	#cd ${outputDir}/anno
	#tRNAscan-SE -qQ -Y -o ${outputDir}/anno/${name}.tRNAscan.tblout -m  ${outputDir}/anno/${name}.tRNAscan.summary -B ${inputFile}

	#perl ${script}/getnum.pl ${inputFile} >${outputDir}/anno/size
	#size=$(cat ${outputDir}/anno/size|awk '{print $1}')
	#/data/public_tools/Infernal/infernal-1.1.2/bin/cmscan -Z ${size} --cut_ga --cpu ${threads} --rfam --nohmmonly --tblout ${outputDir}/anno/${name}.rfam.tblout --fmt 2 --clanin ${database}/Rfam/Rfam.clanin ${database}/Rfam/Rfam.cm ${inputFile}
mkdir -p ${outputDir}/plot
perl ${script}/1.count_genome_GC.pl -i ${inputFile} -o  ${outputDir}/plot/${name} -w ${windowSize} -s ${stepSize}
perl ${script}/2.annotation_formatting.pl  -o ${outputDir}/plot/${name} -prodigalFile ${annopath}/${name}_prodigal/${name}.gff -rnammerFile ${annopath}/${name}_RNAmmer/${name}.RNAmmer.gff2 -trnascanFile ${annopath}/${name}_tRNAscan/${name}.tRNAscan.a -rfamFile /lustre/home/liudongmei/pipeline/pipeline1/database_update/20230625/Rfam/Rfam_entry.tsv -cogPepID2CogID /lustre/home/liudongmei/pipeline/pipeline1/database_update/COG/cog-20.cog.csv -cogCogID2ClassID /lustre/home/liudongmei/pipeline/pipeline1/database_update/COG/cog-20.def.tab -ClassID2FirstID /lustre/home/liudongmei/pipeline/pipeline1/database_update/COG/cog_classify_matched.list -cogFile ${annopath}/${name}_diamond/COG_diamond.txt
perl ${script}/merge.pl ${name} ${outputDir}/plot/ ${inputFile}
perl ${script}/1.count_genome_GC.pl -i ${outputDir}/plot/merge/${name}.fasta -o ${outputDir}/plot/merge/${name} -w ${windowSize} -s ${stepSize}
#perl ${script}/2.annotation_formatting.pl  -o ${outputDir}/plot/${name} -prodigalFile ${outputDir}/anno/${name}.prodigal.gff3 -rnammerFile ${outputDir}/anno/${name}.RNAmmer.gff -trnascanFile ${outputDir}/anno/${name}.tRNAscan.tblout -rfamFile ${outputDir}/anno/${name}.rfam.tblout -rfamRf ${database}/Rfam//Rfam_entry.tsv -cogPepID2CogID ${database}/COG/cog-20.cog.csv -cogCogID2ClassID ${database}/COG/cog-20.def.tab -ClassID2FirstID ${database}/COG/cog_classify_matched.list -cogFile ${outputDir}/anno/${name}.out 

cp ${script}/CGView.config ${outputDir}/plot/merge
sed -i "s#title=#title=\"${name}\"#g" ${outputDir}/plot/merge/CGView.config
contigName=$(awk -F '\t' 'NR==2{print $1}' ${outputDir}/plot/merge/${name}.genomeLength.tab)
mkdir -p  ${outputDir}/plot/merge/CGview/
perl ${script}/3.createCGViewXml.pl -i ${outputDir}/plot/merge/${name} -o  ${outputDir}/plot/merge/${name}.CGView.xml -CGViewFile ${outputDir}/plot/merge/CGView.config -COGColorFile ${script}/COG.color -contigName ${contigName}
java -jar -Djava.awt.headless=true ${script}/cgview.jar -i ${outputDir}/plot/merge/${name}.CGView.xml -s ${outputDir}/plot/merge/CGview/ -R true -e T
