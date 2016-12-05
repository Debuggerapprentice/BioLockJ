

package gMapping;

import bioLockJ.BioLockJExecutor;
import bioLockJ.BioLockJUtils;
import utils.ConfigReader;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;


public class BAM2BCF extends BioLockJExecutor{
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
        BioLockJUtils.requireString(cReader, ConfigReader.BAM_SORTED_DIRECTORY);
        BioLockJUtils.requireString(cReader, ConfigReader.CLUSTER_BATCH_COMMAND);
        BioLockJUtils.requireString(cReader, ConfigReader.REF_SEQUENCE);
        BioLockJUtils.requireString(cReader, ConfigReader.SAMTOOLS_DIRECTORY);
        BioLockJUtils.requireString(cReader, ConfigReader.SCRIPTS_DIRECTORY);
        BioLockJUtils.requireString(cReader, ConfigReader.BCF_DIRECTORY);
        BioLockJUtils.requireString(cReader, ConfigReader.VCF_DIRECTORY);
        BioLockJUtils.requireString(cReader, ConfigReader.BCFTOOLS_DIRECTORY);
    }

    @Override
    public void executeProjectFile(ConfigReader cReader, BufferedWriter logWriter) throws Exception
    {
        this.scripts = new ArrayList<File>();
        File bamSortedDir = BioLockJUtils.requireExistingDirectory(cReader, ConfigReader.BAM_SORTED_DIRECTORY);
        String clusterBatchCommand = BioLockJUtils.requireString(cReader, ConfigReader.CLUSTER_BATCH_COMMAND);
        File refSeq = BioLockJUtils.requireExistingFile(cReader, ConfigReader.REF_SEQUENCE);
        String samtoolsPath = BioLockJUtils.requireString(cReader, ConfigReader.SAMTOOLS_DIRECTORY);
        File scriptDir = BioLockJUtils.requireExistingDirectory(cReader, ConfigReader.SCRIPTS_DIRECTORY);
        File bcfDir = BioLockJUtils.requireExistingFile(cReader, ConfigReader.BCF_DIRECTORY);
        File vcfDir = BioLockJUtils.requireExistingFile(cReader, ConfigReader.VCF_DIRECTORY);
        String bcftoolsPath = BioLockJUtils.requireString(cReader, ConfigReader.BCFTOOLS_DIRECTORY);

        int index =1;
        this.runAllFile = new File(scriptDir.getAbsolutePath() + File.separator + "runAll_" +
                System.currentTimeMillis() + 	".sh");

        BufferedWriter allWriter = new BufferedWriter(new FileWriter(this.runAllFile));

        String[] filesToRun= bamSortedDir.list();

        for( String s : filesToRun)
        {
            File bamFile = new File(bamSortedDir.getAbsolutePath() + File.separator + s);

            if( ! bamFile.isDirectory())
            {
                File script = new File(
                        scriptDir.getAbsolutePath() + File.separator + "run_" + index + "_" +
                                System.currentTimeMillis() + 	"_.sh");

                BufferedWriter writer = new BufferedWriter(new FileWriter(script));

                File outFilebcf = new File( bcfDir+ File.separator + bamFile.getName().replace("_1.fastq_sorted.bam","") + ".bcf");
                File outFilevcf = new File( vcfDir+ File.separator + bamFile.getName().replace("_1.fastq_sorted.bam","") + ".vcf");

                writer.write(samtoolsPath + "samtools mpileup -uf " + refSeq +" "+ bamFile + " | " +
                       bcftoolsPath + "bcftools call --ploidy 1 -mv > " + outFilebcf +"\n");
                writer.write(bcftoolsPath + "bcftools view " + outFilebcf + " > " + outFilevcf);

                File touchFile = new File(script.getAbsolutePath() + FINISHED_SUFFIX );

                if( touchFile.exists())
                    touchFile.delete();

                //writer.write("touch " + touchFile.getAbsolutePath() + "\n");

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

