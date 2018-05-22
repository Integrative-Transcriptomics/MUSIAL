/**
 * 
 */
package main;

import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

import tools.Analyzer;


/**
 * @author Alexander Seitz
 *
 */
public class Musial {

	public static String CLASS_NAME = "title";
	public static String VERSION = "x.x";

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		loadMetadata();
		new Analyzer(args);
			System.out.println("finished");
			System.exit(0);
	}

	private static void loadMetadata(){
		Properties properties = new Properties();
		try {
			//load version
			InputStream in = Musial.class.getResourceAsStream("/version.properties");
			properties.load(in);
			Musial.VERSION = properties.getProperty("version");
			in.close();
			// load title
			in = Musial.class.getResourceAsStream("/title.properties");
			properties.load(in);
			Musial.CLASS_NAME = properties.getProperty("title");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
