package SRATools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.List;

import bioLockJ.BioLockJExecutor;
import bioLockJ.BioLockJUtils;
import utils.ConfigReader;

/*
 * for now, just uses the default kmer size of 31
 */
public class FASTQDumpFromSRRList extends BioLockJExecutor
{
    private File runAllFile = null;
    private List<File> scripts = null;

    @Override
    public File getRunAllFile()
    {
        return runAllFile;
    }

    @Override
    public List<File> getScriptFiles()
    {
        return scripts;
    }

    @Override
    public void checkDependencies(ConfigReader cReader) throws Exception
    {
        BioLockJUtils.requireString(cReader, ConfigReader.FASTQ_DIRECTORY);
        BioLockJUtils.requireString(cReader, ConfigReader.CLUSTER_BATCH_COMMAND);
        BioLockJUtils.requireString(cReader, ConfigReader.BOWTIE2_BINARY);
        BioLockJUtils.requireString(cReader, ConfigReader.SAMTOOLS_BINARY);
        BioLockJUtils.requireString(cReader, ConfigReader.SCRIPTS_DIRECTORY);
        BioLockJUtils.requireString(cReader, ConfigReader.BAM_SORTED_DIRECTORY);
        BioLockJUtils.requireString(cReader, ConfigReader.REF_GENOME_DIRECTORY);
    }

    @Override
    public void executeProjectFile(ConfigReader cReader, BufferedWriter logWriter) throws Exception
    {
        this.scripts = new ArrayList<File>();
        String sraBinaryPath = BioLockJUtils.requireString(cReader, ConfigReader.SRA_BINARY_DIRECTORY);
        String bowtie2Path = BioLockJUtils.requireString(cReader, ConfigReader.BOWTIE2_BINARY);
        String samtoolsPath = BioLockJUtils.requireString(cReader, ConfigReader.SAMTOOLS_BINARY);
        String clusterBatchCommand = BioLockJUtils.requireString(cReader, ConfigReader.CLUSTER_BATCH_COMMAND);
        File bamSortedDir = BioLockJUtils.requireExistingDirectory(cReader, ConfigReader.BAM_SORTED_DIRECTORY);
        File scriptDir = BioLockJUtils.requireExistingDirectory(cReader, ConfigReader.SCRIPTS_DIRECTORY);
        File fastqDir = BioLockJUtils.requireExistingFile(cReader, ConfigReader.FASTQ_DIRECTORY);
        String refDir = BioLockJUtils.requireString(cReader, ConfigReader.REF_GENOME_DIRECTORY);

        int index =1;
        this.runAllFile = new File(scriptDir.getAbsolutePath() + File.separator + "runAll_" +
                System.currentTimeMillis() + 	".sh");

        BufferedWriter allWriter = new BufferedWriter(new FileWriter(this.runAllFile));

        String[] filesToRun= fastqDir.list(new FilenameFilter(){
        public boolean accept(File fastqDir, String fileName) {
                return fileName.endsWith("_1.fastq");
            }
        });

        for( String s : filesToRun)
        {
            File fastqFile = new File(fastqDir.getAbsolutePath() + File.separator + s);

            if( ! fastqFile.isDirectory())
            {
                File script = new File(
                        scriptDir.getAbsolutePath() + File.separator + "run_" + index + "_" +
                                System.currentTimeMillis() + 	"_.sh");

                BufferedWriter writer = new BufferedWriter(new FileWriter(script));

                File outFile = new File( bamSortedDir+ File.separator + fastqFile.getName() + "_sorted.bam");

                writer.write(bowtie2Path + "bowtie2 -x " + refDir + " -1 " +
                        fastqFile.getAbsolutePath() + " -2 " + fastqFile.getAbsolutePath().replace("_1.fastq","_2.fastq") +
                        " | " + samtoolsPath + "samtools view -bS - | " + samtoolsPath + "samtools sort - " + outFile.getAbsolutePath()
                );

                File touchFile = new File(script.getAbsolutePath() + FINISHED_SUFFIX );

                if( touchFile.exists())
                    touchFile.delete();

                writer.write("touch " + touchFile.getAbsolutePath() + "\n");

                writer.flush();  writer.close();
                this.scripts.add(script);

                allWriter.write(clusterBatchCommand + " " + script.getAbsolutePath() + "\n");
                allWriter.flush();

                index++;
            }
        }

        allWriter.flush();  allWriter.close();
    }


}
