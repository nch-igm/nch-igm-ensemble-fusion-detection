SELECT gene1 || '>>' || gene2 AS gene_pair, COUNT(*) AS count
  FROM (SELECT * FROM fusion LEFT JOIN analysis USING (analysis_id)) AS t1
  WHERE tool_id = 2
  GROUP BY gene_pair
  ORDER BY count DESC;
