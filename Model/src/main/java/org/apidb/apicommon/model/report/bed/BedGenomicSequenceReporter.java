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
import org.gusdb.wdk.model.record.attribute.AttributeValue;

public class BedGenomicSequenceReporter extends BedReporter {

  protected List<List<String>> recordAsBedFields(JSONObject config, RecordInstance record){
    return Arrays.asList(recordAsBedFieldsGenomicSequence(config, record));
  }

  private static String stringValue(RecordInstance record, String key){
    try {
      return record.getAttributeValue(key).toString();
    } catch (WdkModelException | WdkUserException e){
      throw new WdkRuntimeException(e);
    }
  }
  
  private static String getSourceId(RecordInstance record){
    return record.getPrimaryKey().getValues().get("source_id");
  }

  private static List<String> recordAsBedFieldsGenomicSequence(JSONObject config, RecordInstance record){
    String featureId = getSourceId(record);
    String chrom = featureId;
    String strand = ".";

    Integer featureLength = Integer.valueOf(stringValue(record, "formatted_length").replaceAll(",", ""));
    Integer segmentStart = config.getInt("start");
    Integer segmentEnd = config.getInt("end");
    String segmentEndDesc = segmentEnd.toString();
    if(segmentEnd == 0){
       segmentEnd = featureLength;
       segmentEndDesc = featureLength + " (end)";
    }
    // Pf3D7_04_v3  | Plasmodium falciparum 3D7 | 1 to 1200490
    StringBuilder sb = new StringBuilder(featureId);
    sb.append("  | ");
    sb.append(stringValue(record, "organism"));
    sb.append(" | (");
    sb.append(""+segmentStart);
    sb.append(" to ");
    sb.append(segmentEndDesc);
    sb.append(")");
    String defline = sb.toString();
    return Arrays.asList(chrom, "" + segmentStart, "" + segmentEnd, defline, ".", strand);
  }

}
