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
public class Musial {//extends Application {

	public static String CLASS_NAME = "title";
	public static String VERSION = "x.x";

	// variables for GUI
//	private Stage primaryStage;
//	private BorderPane rootLayout;
//	private AnchorPane displayLayout;

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		loadMetadata();
		new Analyzer(args);
//		if(args.length<0){
////			launch();
//		}else{
//			if(!(args.length>0)){
//				printHelp();
//				System.exit(1);
//			}
//			String command = args[0];
//			String[] newCommands = Arrays.copyOfRange(args, 1, args.length);
//			Tools toolName = Tools.mvcf2;
//			try{
//				toolName = Tools.valueOf(command);
//			}catch (Exception e){
//				System.err.println("Command not recognized: "+command);
//				printHelp();
//				System.exit(1);
//			}
//			try {
//				Class<?> cl = Class.forName(toolName.getConstructorName());
//				Constructor<?> cons = cl.getConstructor(String[].class);
//				cons.newInstance(new Object[]{newCommands});
//			} catch (ClassNotFoundException | NoSuchMethodException | SecurityException | InstantiationException | IllegalAccessException | IllegalArgumentException | InvocationTargetException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//			try{
//				switch(toolName){
//				case geneAlignment:
//					new GeneAlignment(newCommands);
//					break;
//				case mvcf:
//					new MVCFAnalyzer(newCommands);
//					break;
//				default:
//					printHelp();
//					System.exit(1);
//				}
//			}catch(Exception e){
//				System.err.println(e.getMessage());
//				System.exit(1);
//			}
			System.out.println("finished");
			System.exit(0);
//		}
	}

//	@Override
//	public void start(Stage primaryStage){
//		this.primaryStage = primaryStage;
//		this.primaryStage.setTitle(EAGERTools.CLASS_NAME);
//		initRootLayout();
//		showDisplayLayout();
//	}

//	/**
//	 * 
//	 */
//	private void showDisplayLayout() {
//		try {
//			FXMLLoader loader = new FXMLLoader();
//			loader.setLocation(EAGERTools.class.getResource("/Display.fxml"));
//			this.displayLayout = (AnchorPane) loader.load();
//			DisplayController controller = loader.getController();
//			controller.getGeneAlignmentButton().setTooltip(new Tooltip("Calculate the alignment for certain genes"));
//			controller.getMvcfButton().setTooltip(new Tooltip("run the MultiVCFAnalyzer"));
//			rootLayout.setCenter(this.displayLayout);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//	}

//	/**
//	 * 
//	 */
//	private void initRootLayout() {
//		try {
//			final String os = System.getProperty("os.name");
//			FXMLLoader loader = new FXMLLoader();
//			loader.setLocation(EAGERTools.class.getResource("/RootLayout.fxml"));
//			this.rootLayout = (BorderPane) loader.load();
//			// show the scene containing the root layout
//			Scene scene = new Scene(this.rootLayout);
//			primaryStage.setScene(scene);
//			RootLayoutController controller = loader.getController();
//			if(os != null && os.startsWith("Mac")){
//				controller.getMenuBar().useSystemMenuBarProperty().set(true);
//			}
//			// TODO recent
//			this.primaryStage.show();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//
//	}

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

//	private static void printHelp(){
//		StringBuffer result = new StringBuffer();
//		result.append("\n");
//		result.append("Program: ");
//		result.append(CLASS_NAME);
//		result.append(" (");
//		result.append(DESC);
//		result.append(")");
//		result.append("\n");
//		result.append("Version: ");
//		result.append(VERSION);
//		result.append("\n\n");
//		result.append("Usage:\t");
//		result.append(CLASS_NAME);
//		result.append(" <command> [options]\n\n");
//		result.append("Commands:\n");
//		for(Tools t: Tools.values()){
//			result.append(" ");
//			result.append(t.name());
//			result.append("\t");
//			result.append(t.toString());
//			result.append("\n");
//		}
//		System.err.println(result.toString());
//	}

}
