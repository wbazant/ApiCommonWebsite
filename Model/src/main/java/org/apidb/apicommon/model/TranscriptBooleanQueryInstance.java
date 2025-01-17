package org.apidb.apicommon.model;

import static org.gusdb.fgputil.FormatUtil.NL;

import org.gusdb.fgputil.validation.ValidObjectFactory.RunnableObj;
import org.gusdb.wdk.model.WdkModelException;
import org.gusdb.wdk.model.query.BooleanQueryInstance;
import org.gusdb.wdk.model.query.spec.QueryInstanceSpec;

public class TranscriptBooleanQueryInstance extends BooleanQueryInstance {

  private GeneBooleanQueryInstance genebqi;

  public TranscriptBooleanQueryInstance(RunnableObj<QueryInstanceSpec> spec) {
    super(spec);
    if (!(spec.get().getQuery().get() instanceof TranscriptBooleanQuery)) {
      throw new IllegalStateException("Spec passed to BooleanQueryInstance does not contain a BooleanQuery");
    }
    genebqi = new GeneBooleanQueryInstance(spec);
  }

  @Override
  public String getUncachedSql() throws WdkModelException {

    String booleanGenesSql = genebqi.getUncachedSql();
    boolean pid = TranscriptUtil.isProjectIdInPks(_wdkModel);

    String sql = 
        "SELECT DISTINCT gene_source_id, source_id, " + p(pid,"project_id, ") + "sum(wdk_weight) as wdk_weight, max(left_match) as left_match, max(right_match) as right_match FROM ( " + NL +
        " -- boolean of genes " + NL +
        "WITH genes as (" + booleanGenesSql + ")" + NL +
        " -- major select " + NL +
        "select gene_source_id, source_id, " + p(pid,"project_id, ") + "wdk_weight, decode(sum(left_match), 1, 'Y', 0, 'N') as " + TranscriptBooleanQuery.LEFT_MATCH_COLUMN + ", decode(sum(right_match), 1, 'Y', 0, 'N') as " + TranscriptBooleanQuery.RIGHT_MATCH_COLUMN + NL +
        "from (" + NL +
        "  select left.gene_source_id, left.source_id, " + p(pid,"left.project_id, ") + "genes.wdk_weight, 1 as left_match, 0 as right_match" + NL +
        "  from genes, " + NL +
        "  (" + getLeftSql() + ") left" + NL +
        "  where left.gene_source_id = genes.gene_source_id" + NL +
        "  UNION" + NL +
        "  select right.gene_source_id, right.source_id, " + p(pid,"right.project_id, ") + "genes.wdk_weight, 0 as left_match, 1 as right_match" + NL +
        "  from genes, " + NL +
        "  (" + getRightSql() + ") right" + NL +
        "  where right.gene_source_id = genes.gene_source_id" + NL +
        "  UNION" + NL +
        "  select ta.gene_source_id, ta.source_id, " + p(pid,"genes.project_id, ") + "genes.wdk_weight, 0 as left_match, 0 as right_match" + NL +
        "  from genes, apidbtuning.transcriptattributes ta" + NL +
        "  where genes.gene_source_id = ta.gene_source_id) big" + NL +
        "group by (gene_source_id, source_id, " + p(pid,"project_id, ") + "wdk_weight)" +
        ") group by (gene_source_id, source_id" + p(pid,", project_id") + ")";
    return sql;

  }

  // returns the passed string if projectIdRequired is true, else empty string
  private String p(boolean projectIdRequired, String string) {
    return projectIdRequired ? string : "";
  }
}
