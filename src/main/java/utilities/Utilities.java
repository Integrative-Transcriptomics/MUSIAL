/**
 * 
 */
package utilities;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.LinkedList;
import java.util.List;

/**
 * @author Alexander Seitz
 *
 */
public class Utilities {
	
	public static void copyFileUsingIO(File sourceFile, File destinationFile) {
		InputStream inputStreamData = null;
		OutputStream outputStreamData = null;

		try {
			inputStreamData = new BufferedInputStream(new FileInputStream(sourceFile));
			outputStreamData = new BufferedOutputStream(new FileOutputStream(destinationFile));
			byte[] buffer = new byte[1024];
			int length;
			while ((length = inputStreamData.read(buffer)) > 0) {
				outputStreamData.write(buffer, 0, length);
			}
		}catch(IOException e) {
			System.err.println("error copying file");
			e.printStackTrace();
				
		} finally {
			try {
				inputStreamData.close();
				outputStreamData.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	public static boolean fileExists(String f){
		return new File(f).exists();
	}
	
	public static List<String> readFileLinewise(String filename){
		List<String> result = new LinkedList<String>();
		try {
			@SuppressWarnings("resource")
			BufferedReader br = new BufferedReader(new FileReader(filename));
			String line = "";
			while((line = br.readLine()) != null){
				result.add(line);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return result;
	}
	
	public static void writeToFile(String s, File f){
		try {
			PrintWriter outMerged = new PrintWriter(new BufferedWriter(new FileWriter(f, false)));
			outMerged.println(s);
			outMerged.flush();
			outMerged.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static String removeTrailingSlashFromFolder(String name) {
		if (name.endsWith("/")) {
			name = name.substring(0, name.length() - 1);
		}
		return name;
	}
	
	public static void createOutFolder(String folder) {
		String[] createOutputFolder = { "mkdir", "-p", folder };
		runCommand(createOutputFolder, "", "", "");
	}
	
	
	public static void runScript(String script) {
		String[] runScript = {"sh", script};
		runCommand(runScript, "", "", "");
	}
	
	public static void runCommand(String[] command, String errorLog, String output, String runInDir){
		try {
			ProcessBuilder pb = new ProcessBuilder(command);
			if(errorLog.length()>0){
				pb = pb.redirectError(new File(errorLog));
			}
			if(output.length()>0){
				pb = pb.redirectOutput(new File(output));
			}
			if(runInDir.length()>0) {
				pb.directory(new File(runInDir));
			}
			Process process = pb.start();
			process.waitFor();
		} catch (IOException | InterruptedException e) {
		}
	}
	
	
	public static void generateSymlink(String file, String link) {
		String[] generateSymlink = { "ln", "-s", file, link };
		runCommand(generateSymlink, "", "", "");
	}
	
	
	public static String getCurrentDateTime() {
		Calendar cal = Calendar.getInstance();
		SimpleDateFormat sdf = new SimpleDateFormat("yyyyMMdd-HHmmss");
		return sdf.format(cal.getTime());
	}
	

	
	public static Boolean isInteger(String s){
		try{
			Integer.parseInt(s);
			return true;
		}catch(Exception e){
			return false;
		}
	}
	
	public static Boolean isDouble(String s){
		try{
			Double.parseDouble(s);
			return true;
		}catch(Exception e){
			return false;
		}
	}
	
	public static List<File> getEagerSamples(File eagerFolder) {
		List<File> eagerSamples = new LinkedList<File>();
		if(!Utilities.fileExists(eagerFolder.getAbsolutePath())){
			System.err.println("Could not fine the Eager input folder");
		}
		eagerSamples = new LinkedList<File>();
		if(isEagerOutputFolder(eagerFolder.getAbsolutePath())){
			eagerSamples.add(eagerFolder);
		}
		String[] possibleSamples = eagerFolder.list();
		for(String possibleSample: possibleSamples){
			String possSamp = eagerFolder.getAbsolutePath() + "/" + possibleSample;
			if(isEagerOutputFolder(possSamp)){
				eagerSamples.add(new  File(possSamp));
			}
		}
		return eagerSamples;
	}

	/**
	 * @param inputFilder
	 * check if given folder is an EAGER output folder
	 * @return true if the given folder is an EAGER output folder, else false
	 */
	public static boolean isEagerOutputFolder(String inputFilder) {
		File folder = new File(inputFilder);
		if(!folder.exists()){
			return false;
		}
		if(!folder.isDirectory()){
			return false;
		}
		String[] names = folder.list();
		for(String name: names){
			File currFile = new File(folder.getAbsolutePath()+"/"+name);
			if(currFile.isDirectory()){
				for(int i=0; i<=12; i++){
					if(currFile.getName().startsWith(""+i+"-")){
						return true;
					}
				}
			}
		}
		return false;
	}
	
	public static void ExportResource(String resourceName, String outputFile) {
        InputStream stream = null;
        OutputStream resStreamOut = null;
        try {
            stream = Utilities.class.getResourceAsStream(resourceName);//note that each / is a directory down in the "jar tree" been the jar the root of the tree
            if(stream == null) {
            	System.err.println("Cannot get resource \"" + resourceName + "\" from Jar file.");
            	return;
            }

            int readBytes;
            byte[] buffer = new byte[4096];
            resStreamOut = new FileOutputStream(outputFile);
            while ((readBytes = stream.read(buffer)) > 0) {
                resStreamOut.write(buffer, 0, readBytes);
            }
        } catch (Exception ex) {
        	ex.printStackTrace();
        	return;
        } finally {
            try {
				stream.close();
				resStreamOut.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
        }
    }

}
