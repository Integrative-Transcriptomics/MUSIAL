package datastructure;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import exceptions.MusialIOException;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ConcurrentSkipListSet;
import main.Musial;

public class VariantsDictionary {

  public final VariantsDictionaryParameters parameters;

  public final HashMap<String, FeatureEntry> features = new HashMap<>();

  public final HashMap<String, SampleEntry> samples = new HashMap<>();

  public final String software = Musial.NAME + Musial.VERSION;

  public final String date = new SimpleDateFormat("dd/MM/yyyy").format(new Date());

  public final String chromosome;

  public final ConcurrentSkipListMap<Integer, ConcurrentSkipListMap<String, NucleotideVariantAnnotationEntry>>
      variants =
      new ConcurrentSkipListMap<>(Integer::compare);

  public transient ConcurrentSkipListSet<String> novelVariants = new ConcurrentSkipListSet<>((s1, s2) -> {
    int p1 = Integer.parseInt(s1.split("@")[0]);
    int p2 = Integer.parseInt(s2.split("@")[0]);
    if (p1 != p2) {
      return Integer.compare(p1, p2);
    } else {
      String c1 = s1.split("@")[2];
      String c2 = s2.split("@")[2];
      return c1.compareTo(c2);
    }
  });

  public VariantsDictionary(Double minCoverage, Double minFrequency, Double minHet, Double maxHet, Double minQuality,
                            String chromosome) {
    this.parameters = new VariantsDictionaryParameters(minCoverage, minFrequency, minHet, maxHet, minQuality);
    this.chromosome = chromosome;
  }

  public void addVariant(String featureId, int referencePosition, String variantContent,
                         String referenceContent, String sampleId,
                         boolean isPrimary, boolean isRejected, double quality, double coverage, double frequency) {
    // Update variant information in `this.variants`.
    if (!this.variants.containsKey(referencePosition)) {
      this.variants.put(referencePosition, new ConcurrentSkipListMap<>());
    }
    if (!this.variants.get(referencePosition).containsKey(variantContent)) {
      this.variants.get(referencePosition).put(variantContent, new NucleotideVariantAnnotationEntry());
      this.novelVariants.add(referencePosition + "@" + referenceContent + "@" + variantContent);
    }
    if (!this.variants.get(referencePosition).get(variantContent).occurrence.containsKey(sampleId)) {
      this.variants.get(referencePosition).get(variantContent).occurrence.put(sampleId,
          NucleotideVariantAnnotationEntry
              .constructSampleSpecificAnnotation(isRejected, isPrimary, quality, frequency, coverage)
      );
      if (isPrimary && !isRejected) {
        if (!samples.get(sampleId).annotations.containsKey(featureId + "vSwab")) {
          this.samples.get(sampleId).annotations.put(featureId + "vSwab", variantContent + "@" + referencePosition);
        } else {
          this.samples.get(sampleId).annotations
              .put(featureId + "vSwab", this.samples.get(sampleId).annotations.get(featureId + "vSwab")
                  .concat("|" + variantContent + "@" + referencePosition)
              );
        }
      }
    }
  }

  public void addVariantAnnotation(int position, String content, Map<String, String> annotations) {
    if (this.variants.containsKey(position) && this.variants.get(position).containsKey(content)) {
      this.variants.get(position).get(content).annotations.putAll(annotations);
    }

  }

  public void dump(File outfile) throws MusialIOException {
    if (!outfile.exists()) {
      throw new MusialIOException("The specified output file does not exist:\t" + outfile.getAbsolutePath());
    }
    try {
      Gson gson;
      String dumpString;
      String dumpFilePath = outfile.getAbsolutePath();

      /* Write compressed JSON.
      gson = new GsonBuilder().create();
      dumpString = gson.toJson(this);
      try (FileOutputStream output =
               new FileOutputStream(dumpFilePath.endsWith(".gz") ? dumpFilePath : dumpFilePath + ".gz");
           Writer writer = new OutputStreamWriter(new GZIPOutputStream(output), StandardCharsets.UTF_8)) {
        writer.write(dumpString);
      }*/

      // Write to pretty JSON.
      gson = new GsonBuilder().setPrettyPrinting().create();
      dumpString = gson.toJson(this);
      FileWriter variantDBWriter = new FileWriter(dumpFilePath.endsWith(".gz") ? dumpFilePath.replace(".gz", "") :
          dumpFilePath);
      variantDBWriter.write(dumpString);
      variantDBWriter.close();
    } catch (IOException e) {
      throw new MusialIOException("Failed to write to output file:\t" + outfile.getAbsolutePath());
    }
  }

  public String getNucleotideSequence(String fId, String sId) {
    if (sId.equals("REFERENCE")) {
      return features.get(fId).nucleotideSequence;
    }
    if (!samples.get(sId).annotations.containsKey(fId + "vSwab")) {
      return null;
    }
    char[] referenceSequence = features.get(fId).nucleotideSequence.toCharArray();
    StringBuilder sampleSequenceBuilder = new StringBuilder();
    HashMap<Integer, String> sampleVariants = new HashMap<>();
    String variantContent;
    int referencePosition;
    int referenceStart;
    int d;
    for (String vInformation : samples.get(sId).annotations.get(fId + "vSwab").split("\\|")) {
      sampleVariants.put(Integer.parseInt(vInformation.split("@")[1]), vInformation.split("@")[0]);
    }
    referenceStart = features.get(fId).start;
    d = 0; // Counter for deleted positions.
    for (int i = 0; i < referenceSequence.length; i++) {
      referencePosition = referenceStart + i;
      if (d > 0) {
        d--;
      }
      if (sampleVariants.containsKey(referencePosition)) {
        variantContent = sampleVariants.get(referencePosition);
        if (variantContent.contains("-")) {
          // Variant is a deletion.
          d = variantContent.length();
          sampleSequenceBuilder.append(referenceSequence[i]);
        } else {
          sampleSequenceBuilder.append(variantContent);
        }
      } else if (d == 0) {
        sampleSequenceBuilder.append(referenceSequence[i]);
      }
    }
    return sampleSequenceBuilder.toString();
  }

  public HashMap<String, ArrayList<String>> getNucleotideSequences(String fId) {
    FeatureEntry featureEntry = features.get(fId);
    HashMap<String, ArrayList<String>> nucleotideSequenceMap = new HashMap<>();
    ArrayList<String> samplesIndexList = new ArrayList<>();
    String sampleSequence;
    nucleotideSequenceMap.put(
        featureEntry.nucleotideSequence,
        new ArrayList<>() {{
          add("REFERENCE");
        }}
    );
    for (String sId : samples.keySet()) {
      sampleSequence = getNucleotideSequence(fId, sId);
      if (sampleSequence == null) {
        continue;
      }
      if (nucleotideSequenceMap.containsKey(sampleSequence)) {
        nucleotideSequenceMap.get(sampleSequence).add(sId);
      } else {
        samplesIndexList.add(sId);
        nucleotideSequenceMap.put(sampleSequence, samplesIndexList);
        samplesIndexList = new ArrayList<>();
      }
    }
    return nucleotideSequenceMap;
  }
}
