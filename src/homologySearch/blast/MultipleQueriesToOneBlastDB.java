package homologySearch.blast;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import bioLockJ.BioJLockUtils;
import bioLockJ.BioLockJExecutor;
import homologySearch.BreakUpFastaSequence;
import utils.ConfigReader;

/**
 * Takes in a BLAST_QUERY_DIRECTORY that should only contain FASTA files (subdirectories are ignored)
 * Takes in a FASTA_FILE_TO_FORMAT_FOR_BLAST_DB (assuming formatted for example by FormatSingleBlastDatabase)
 * Writes multiple scripts to SCRIPTS_DIR_FOR_BLAST_QUERY
 * Writes results to BLAST_OUTPUT_DIRECTORY
 * 
 * will issue BLAST_PRELIMINARY_STRING if defined
 * requires BLAST_BIN_DIR to be defined
 * 
 * CLUSTER_BATCH_COMMAND must be defined (e.g. qsub -q "viper" ) where viper is the name of the cluster
 */
public class MultipleQueriesToOneBlastDB extends BioLockJExecutor
{
	private File runAllFile = null;
	private List<File> scripts = null;
	
	@Override
	public void executeProjectFile(File projectFile) throws Exception
	{
		this.scripts = new ArrayList<File>();
		ConfigReader cReader = new ConfigReader(projectFile);
		String blastBinDin = BioJLockUtils.requireString(cReader, ConfigReader.BLAST_BINARY_DIR);
		File blastQueryDir = BioJLockUtils.requireExistingDirectory(cReader, ConfigReader.BLAST_QUERY_DIRECTORY);
		File blastDatabaseFile = BioJLockUtils.requireExistingFile(cReader, ConfigReader.FASTA_FILE_TO_FORMAT_FOR_BLAST_DB);
		File blastOutputDirectory = BioJLockUtils.requireExistingDirectory(cReader, ConfigReader.BLAST_OUTPUT_DIRECTORY);
		
		File scriptDir = BioJLockUtils.requireExistingDirectory(cReader, ConfigReader.SCRIPTS_DIR_FOR_BLAST_QUERY);
	
		File logDir = BioJLockUtils.createLogDirectory(scriptDir, BreakUpFastaSequence.class.getSimpleName());
		BioJLockUtils.copyPropertiesFile(projectFile, logDir);
		
		BufferedWriter logWriter = new BufferedWriter(new FileWriter(new File(
				logDir.getAbsolutePath() + File.separator + BreakUpFastaSequence.class.getSimpleName() 
				 +"log.txt")));
		
		String clusterBatchCommand = BioJLockUtils.requireString(cReader, ConfigReader.CLUSTER_BATCH_COMMAND);
		
		int index =1;
		this.runAllFile = new File(scriptDir.getAbsolutePath() + File.separator + "runAll_" + 
				System.currentTimeMillis() + 	".sh");
		
		BufferedWriter allWriter = new BufferedWriter(new FileWriter(this.runAllFile));
		
		String[] filesToFormat = blastQueryDir.list();
		
		for( String s : filesToFormat)
		{
			File fastaFile = new File(blastQueryDir.getAbsolutePath() + File.separator + s);
			
			if( ! fastaFile.isDirectory())
			{
				File script = new File(
						scriptDir.getAbsolutePath() + File.separator + "run_" + index + "_" +
								System.currentTimeMillis() + 	"_.sh");
				
				BufferedWriter writer = new BufferedWriter(new FileWriter(script));
				
				String prelimString = cReader.getAProperty(ConfigReader.BLAST_PRELIMINARY_STRING);
				
				if( prelimString != null)
					writer.write(prelimString + "\n");
				
				File outFile = new File( blastOutputDirectory + File.separator + fastaFile.getName() + "_to_" + 
							blastDatabaseFile.getName() + ".txt");
				
				writer.write(blastBinDin + "/blastn -db " + 
						blastDatabaseFile.getAbsolutePath() + " -out " + 
							outFile.getAbsolutePath() +  
							" -query " +fastaFile.getAbsolutePath() + 
							" -outfmt 6\n");
				
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

		logWriter.write("successful completion at " + new Date().toString() + "\n"); 
		logWriter.flush(); logWriter.close();
		BioJLockUtils.appendSuccessToPropertyFile(projectFile, this.getClass().getName(), logDir);
	}
	
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
	
	public static void main(String[] args) throws Exception
	{
		if( args.length != 1)
		{
			System.out.println("Usage " + MultipleQueriesToOneBlastDB.class.getName() + " pathToPropertyFile" );
			System.exit(1);
		}
		
		File propFile = BioJLockUtils.findProperyFile(args);
		new MultipleQueriesToOneBlastDB().executeProjectFile(propFile);
	}
}