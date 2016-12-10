package gMapping;

import bioLockJ.BioLockJExecutor;
import bioLockJ.BioLockJUtils;
import utils.ConfigReader;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.List;


public class GzipAllVCF extends BioLockJExecutor
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
        BioLockJUtils.requireString(cReader, ConfigReader.CLUSTER_BATCH_COMMAND);
        BioLockJUtils.requireString(cReader, ConfigReader.SCRIPTS_DIRECTORY);
        BioLockJUtils.requireString(cReader, ConfigReader.VCF_DIRECTORY);
        BioLockJUtils.requireString(cReader, ConfigReader.HTSLIB_DIRECTORY);
        BioLockJUtils.requireString(cReader, ConfigReader.VCF_GZIP_DIRECTORY);
    }

    @Override
    public void executeProjectFile(ConfigReader cReader, BufferedWriter logWriter) throws Exception
    {
        this.scripts = new ArrayList<File>();
        String clusterBatchCommand = BioLockJUtils.requireString(cReader, ConfigReader.CLUSTER_BATCH_COMMAND);
        File scriptDir = BioLockJUtils.requireExistingDirectory(cReader, ConfigReader.SCRIPTS_DIRECTORY);
        File vcfDir = BioLockJUtils.requireExistingDirectory(cReader, ConfigReader.VCF_DIRECTORY);
        String htslibDir = BioLockJUtils.requireString(cReader, ConfigReader.HTSLIB_DIRECTORY);
        File vcfGzipDir = BioLockJUtils.requireExistingDirectory(cReader, ConfigReader.VCF_GZIP_DIRECTORY);

        int index =1;
        this.runAllFile = new File(scriptDir.getAbsolutePath() + File.separator + "runAll_" +
                System.currentTimeMillis() + 	".sh");

        BufferedWriter allWriter = new BufferedWriter(new FileWriter(this.runAllFile));

        String[] filesToRun= vcfDir.list(new FilenameFilter(){
        public boolean accept(File vcfDir, String fileName) {
                return fileName.endsWith(".vcf");
            }
        });

        for( String s : filesToRun)
        {
            File vcfFile = new File(vcfDir.getAbsolutePath() + File.separator + s);

            if( ! vcfFile.isDirectory())
            {
                File script = new File(
                        scriptDir.getAbsolutePath() + File.separator + "run_" + index + "_" +
                                System.currentTimeMillis() + 	"_.sh");

                BufferedWriter writer = new BufferedWriter(new FileWriter(script));

                File outFileGz = new File( vcfGzipDir + File.separator + vcfFile.getName() + ".gz");
                //File outFileTbi = new File( vcfGzipDir + File.separator + vcfFile.getName() + ".gz.tbi";


                writer.write( htslibDir + "bgzip " + vcfFile + " > " + outFileGz);

                writer.write( htslibDir + "tabix -p vcf " + outFileGz);

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
