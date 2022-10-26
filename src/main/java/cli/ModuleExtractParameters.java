package cli;

import com.google.gson.internal.LinkedTreeMap;
import components.Logging;
import components.Validation;
import exceptions.MusialException;

import java.io.File;
import java.util.ArrayList;

/**
 * Parses command line interface arguments for the `MUSIAL EXTRACT*` modules.
 * <p>
 * Used to parse and validate passed command line options or print help information to the user. Parsed arguments are
 * stored to be accessible for other components of the tool.
 *
 * @author Simon Hackl
 * @version 2.1
 * @since 2.1
 */
public class ModuleExtractParameters {
    /**
     * {@link File} object specifying the input file.
     */
    public File inputFile;
    /**
     * {@link File} object specifying the output directory.
     */
    public File outputDirectory;
    /**
     * {@link ArrayList} of samples to consider.
     */
    public ArrayList<String> samples;
    /**
     * {@link ArrayList} of features to consider.
     */
    public ArrayList<String> features;
    /**
     * {@link Boolean} whether to include only SNVs.
     */
    public boolean SNVOnly;

    @SuppressWarnings("FieldCanBeLocal")
    private final String EXCEPTION_PREFIX = "(Module BUILD Configuration)";

    /**
     * Constructor of the {@link ModuleExtractParameters} class.
     * <p>
     * Used to parse and validate command line interface input. Parsed arguments are stored - and accessible by other
     * components - via class properties.
     *
     * @param parameters {@link LinkedTreeMap} containing parameter information.
     * @throws MusialException If IO file validation fails; If CLI parameter validation fails.
     */
    public ModuleExtractParameters(LinkedTreeMap<Object, Object> parameters)
            throws MusialException {
        // Parse input file
        if (!parameters.containsKey("inputFile")) {
            throw new MusialException(EXCEPTION_PREFIX + " Missing " + Logging.colorParameter("inputFile") + "; expected path to file.");
        }
        File i = new File((String) parameters.get("inputFile"));
        if (Validation.isDirectory(new File(i.getParent()))) {
            this.inputFile = i;
        } else {
            throw new MusialException(EXCEPTION_PREFIX + " Invalid " + Logging.colorParameter("inputFile") + " " + Logging.colorParameter((String) parameters.get("input")) + "; unable to access directory.");
        }

        // Parse output directory
        if (!parameters.containsKey("outputDirectory")) {
            throw new MusialException(EXCEPTION_PREFIX + " Missing " + Logging.colorParameter("outputDirectory") + "; expected path to directory.");
        }
        File o = new File((String) parameters.get("outputDirectory"));
        if (Validation.isDirectory(o)) {
            this.outputDirectory = o;
        } else {
            throw new MusialException(EXCEPTION_PREFIX + " Invalid " + Logging.colorParameter("outputDirectory") + " " + Logging.colorParameter((String) parameters.get("outputDirectory")) + ".");
        }

        // Parse samples to consider
        if (parameters.containsKey("samples")) {
            //noinspection unchecked
            this.samples = (ArrayList<String>) parameters.get("samples");
        } else {
            this.samples = new ArrayList<>();
        }

        // Parse features to consider
        if (parameters.containsKey("features")) {
            //noinspection unchecked
            this.features = (ArrayList<String>) parameters.get("features");
        } else {
            this.features = new ArrayList<>();
        }

        // Parse boolean parameter whether to exclude indels.
        if (parameters.containsKey("snvOnly")) {
            this.SNVOnly = (boolean) parameters.get("snvOnly");
        } else {
            this.SNVOnly = false;
        }
    }

}
