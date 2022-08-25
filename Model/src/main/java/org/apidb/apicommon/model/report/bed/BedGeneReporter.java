package org.gusdb.wdk.model.report.reporter.bed;

import org.gusdb.wdk.model.report.reporter.bed.BedReporter;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import org.gusdb.wdk.model.record.RecordInstance;
import org.gusdb.wdk.model.record.TableValue;
import org.json.JSONObject;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.gusdb.wdk.model.WdkModelException;
import org.gusdb.wdk.model.WdkUserException;
import org.gusdb.wdk.model.WdkRuntimeException;
import org.gusdb.wdk.model.report.ReporterConfigException;
import org.gusdb.wdk.model.record.attribute.AttributeValue;

import org.apidb.apicommon.model.TranscriptUtil;

public class BedGeneReporter extends BedReporter {

  private String _originalQuestionName;

  @Override
  public BedGeneReporter configure(JSONObject config) throws ReporterConfigException {
    try {
      _originalQuestionName = _baseAnswer.getAnswerSpec().getQuestion().getName();
      if(configNeedsGeneAnswer(config)){
        _baseAnswer = TranscriptUtil.transformToGeneAnswer(_baseAnswer);
      }
      super.configure(config);
      return this;
    }
    catch (WdkUserException e) {
      throw new ReporterConfigException(e.getMessage());
    }
    catch (WdkModelException e) {
      throw new WdkRuntimeException("Could not create in-memory step from incoming answer spec", e);
    }
  }

  @Override
  public String getDownloadFileName() {
    return _originalQuestionName + ".tsv";
  }

  private static boolean configNeedsGeneAnswer(JSONObject config){
    String type = config.getString("type");
    switch(type){
      case "genomic":
      case "protein":
        return false;
      case "protein_features":
      case "genomic_features":
        return true;
      default:
        throw new WdkRuntimeException(String.format("Unknown sequence type: %s", type));
    }
  }

  protected List<List<String>> recordAsBedFields(JSONObject config, RecordInstance record){
    String type = config.getString("type");
    switch(type){
      case "genomic":
        return Arrays.asList(recordAsBedFieldsGenomic(config, record));
      case "protein":
        return Arrays.asList(recordAsBedFieldsProtein(config, record));
      case "protein_features":
        String proteinFeature = config.getString("proteinFeature");
        switch (proteinFeature){
          case "interpro":
            return featuresAsListOfBedFieldsProteinInterpro(config, record);
          case "signalp":
          case "tmhmm":
            throw new WdkRuntimeException(String.format("Unknown protein feature: %s", proteinFeature));
        }
      case "genomic_features":
        String genomicFeature = config.getString("genomicFeature");
        switch (genomicFeature){
          case "low_complexity_regions":
            throw new WdkRuntimeException("TODO low_complexity_regions");
          default:
            throw new WdkRuntimeException(String.format("Unknown genomic feature: %s", genomicFeature));
        }
      default:
        throw new WdkRuntimeException(String.format("Unknown sequence type: %s", type));
    }
  }


  // "PbANKA_01_v3:438265..440094(-)"
  private static Pattern locationTextPattern = Pattern.compile("^(.*):(\\d+)..(\\d+)\\((\\+|-)\\)$");

  private static String stringValue(RecordInstance record, String key){
    try {
      return record.getAttributeValue(key).toString();
    } catch (WdkModelException | WdkUserException e){
      throw new WdkRuntimeException(e);
    }
  }

  private static String shortSign(String sign){
    switch(sign){
      case "plus":
        return "+";
      case "minus":
        return "-";
      default:
        return sign;
    }
  }
  private static String longStrand(String shortSign){
    switch(shortSign){
      case "+":
        return "forward";
      case "-":
        return "reverse";
      default:
        return shortSign;
    }
  }

  private static Integer getPositionGenomic(JSONObject config, Integer featureStart, Integer featureEnd, 
      String offsetKey, String signKey, String anchorKey){
    String sign = config.getString(signKey);
    Integer offset;
    switch(sign){
      case "plus":
        offset = config.getInt(offsetKey);
        break;
      case "minus":
        offset = - config.getInt(offsetKey);
        break;
      default:
        throw new WdkRuntimeException(String.format("%s value should be 'plus' or 'minus', got: %s", signKey, sign));
    }

    Integer result;
    String anchor = config.getString(anchorKey);
    switch(anchor){
      case "Start":
        result = featureStart + offset;
        break;
      case "End":
        result = featureEnd + offset;
        break;
      default:
        throw new WdkRuntimeException(String.format("%s value should be 'Start' or 'End', got: %s", anchorKey, anchor));
    }
    return result;
  }

  private static Integer getPositionProtein(JSONObject config, Integer featureLength, 
      String offsetKey, String anchorKey){
    Integer offset = config.getInt(offsetKey);

    Integer result;
    String anchor = config.getString(anchorKey);
    switch(anchor){
      case "Start":
        result = offset;
        break;
      case "End":
        result = featureLength - offset;
        break;
      default:
        throw new WdkRuntimeException(String.format("%s value should be 'Start' or 'End', got: %s", anchorKey, anchor));
    }
    return result;
  }

  private static String getPositionDescGenomic(JSONObject config,
      String offsetKey, String signKey, String anchorKey){
    Integer offset = Integer.valueOf(config.getInt(offsetKey));
    if(offset == 0){
        return config.getString(anchorKey);
    } else {
        return config.getString(anchorKey) + shortSign(config.getString(signKey)) + offset.toString();
    }
  }
  private static String getPositionDescProtein(JSONObject config,
      String offsetKey, String sign, String anchorKey){
    Integer offset = Integer.valueOf(config.getInt(offsetKey));
    if(offset == 0){
        return config.getString(anchorKey);
    } else {
        return config.getString(anchorKey) + sign + offset.toString();
    }
  }

  private static Matcher matchLocationCoords(RecordInstance record, String key, Pattern p){
    String text = stringValue(record, key);
    Matcher m = p.matcher(text);
    if (!m.matches()){
      throw new WdkRuntimeException(String.format("attribute %s with value %s not matching pattern %s", key, text, p.toString()));
    }
    return m;
  }
  
  private static String getSourceId(RecordInstance record){
    return record.getPrimaryKey().getValues().get("source_id");
  }

  private static List<String> recordAsBedFieldsGenomic(JSONObject config, RecordInstance record){
    String featureId = getSourceId(record);
    Matcher m = matchLocationCoords(record, "location_text", locationTextPattern);
    String chrom = m.group(1);
    Integer featureStart = Integer.valueOf(m.group(2));
    Integer featureEnd = Integer.valueOf(m.group(3));
    String strand = m.group(4);
    String chromDesc = "genomic | %s %s".format(chrom, longStrand(strand));
    String featureType = "gene";
    Integer segmentStart = getPositionGenomic(config, featureStart, featureEnd, "upstreamOffset", "upstreamSign", "upstreamAnchor");
    Integer segmentEnd = getPositionGenomic(config, featureStart, featureEnd, "downstreamOffset", "downstreamSign", "downstreamAnchor");
    
    String defline;
    StringBuilder sb = new StringBuilder(featureId);
    if("1".equals(config.getString("onlyIdDefLine"))){
      // PBANKA_0111300
      defline = sb.toString();
    } else {
      /*
       * was:
       PBANKA_0111300  | Plasmodium berghei ANKA | conserved Plasmodium protein, unknown function | genomic | PbANKA_01_v3 reverse | (geneStart+0 to geneEnd+0)| length=1830
       aiming to change:
       - omit "+0" or "-0" in feature description
       - extra space after region description
       */
      sb.append("  | ");
      sb.append(stringValue(record, "organism"));
      sb.append(" | ");
      sb.append(stringValue(record, "gene_product"));
      sb.append(" | ");
      sb.append(chromDesc);
      sb.append(" | (");
      sb.append(featureType);
      sb.append(getPositionDescGenomic(config, "upstreamOffset", "upstreamSign", "upstreamAnchor"));
      sb.append(" to ");
      sb.append(featureType);
      sb.append(getPositionDescGenomic(config, "downstreamOffset", "downstreamSign", "downstreamAnchor"));
      sb.append(")| length=");
      sb.append(""+(segmentEnd - segmentStart));
      defline = sb.toString();
    }
    return Arrays.asList(chrom, "" + segmentStart, "" + segmentEnd, defline, ".", strand);
  }

  private static List<String> recordAsBedFieldsProtein(JSONObject config, RecordInstance record){
    String featureId = getSourceId(record);
    String chrom = featureId;
    Integer featureLength = Integer.valueOf(stringValue(record, "protein_length"));
    String strand = ".";
    String chromDesc = "protein";
    String featureType = "protein";

    Integer segmentStart = getPositionProtein(config, featureLength, "startOffset3", "startAnchor3");
    Integer segmentEnd = getPositionProtein(config, featureLength, "endOffset3", "endAnchor3");
    
    String defline;
    StringBuilder sb = new StringBuilder(featureId);
    if("1".equals(config.getString("onlyIdDefLine"))){
      // PBANKA_0111300
      defline = sb.toString();
    } else {
      /*
       * was:
       PBANKA_0111300  | Plasmodium berghei ANKA | conserved Plasmodium protein, unknown function | protein | length=85
       aiming to change:
       - add protein region descriptions
       */
      sb.append("  | ");
      sb.append(stringValue(record, "organism"));
      sb.append(" | ");
      sb.append(stringValue(record, "gene_product"));
      sb.append(" | ");
      sb.append(chromDesc);
      sb.append(" | (");
      sb.append(featureType);
      sb.append(getPositionDescProtein(config, "startOffset3", "+", "startAnchor3"));
      // #TODO descs
      sb.append(" to ");
      sb.append(getPositionDescProtein(config, "endOffset3", "-", "endAnchor3"));
      sb.append(")| length=");
      sb.append(""+(segmentEnd - segmentStart));
      defline = sb.toString();
    }
    return Arrays.asList(chrom, "" + segmentStart, "" + segmentEnd, defline, ".", strand);
  }

  private static List<List<String>> featuresAsListOfBedFieldsProteinInterpro(JSONObject config, RecordInstance record){
    String featureId = getSourceId(record);
    List<List<String>> result = new ArrayList<>();
    try {
      TableValue interproRows = record.getTableValue("InterPro");
      for (Map<String, AttributeValue> interproRow : interproRows) {
        String start = interproRow.get("interpro_start_min").toString();
        String end = interproRow.get("interpro_end_min").toString();
        /*
         * The start and end coordinates are on the protein,
         * but for identical positions, the rows get merged.
         * Hence 'transcript_ids' that we split.
         * Example: PF3D7_0108400
         */
        for(String transcriptId: interproRow.get("transcript_ids").toString().split(", ")){
          String chrom = transcriptId;
          String defline = transcriptId + "::" + interproRow.get("interpro_primary_id").toString();
          result.add(Arrays.asList(chrom, start, end, defline, ".", "."));
        }
      }
    } catch (WdkModelException | WdkUserException e){
      throw new WdkRuntimeException(e);
    }

    return result;
  }
}
