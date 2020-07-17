create table tool (
  tool_id serial primary key,
  tool_name varchar(255) not null
);

create table subject (
  subject_id serial primary key,
  subject_name varchar(255) not null,
  phenotype varchar(2048)
);
alter table subject add constraint subject_name_unique unique (subject_name);

create table sample (
  sample_id serial primary key,
  subject_id integer not null references subject,
  sample_name varchar(255) not null,
  storage_method varchar(255),
  num_reads integer
);
alter table sample add constraint sample_name_unique unique (sample_name);

create table analysis (
  analysis_id serial primary key,
  sample_id integer not null references sample,
  tool_id integer not null references tool,
  tool_version varchar(255),
  unique (sample_id, tool_id)
);

create table fusion (
  fusion_id serial primary key,
  analysis_id integer not null references analysis,
  fusion_rank integer not null,
  gene1 varchar(2048) not null,
  chrom1 varchar(255) not null,
  break1 integer,
  strand1_plus boolean,
  gene2 varchar(2048) not null,
  chrom2 varchar(255) not null,
  break2 integer,
  strand2_plus boolean
);
alter table fusion add constraint analysis_rank_unique unique (analysis_id, fusion_rank);
ALTER TABLE fusion ALTER COLUMN strand1_plus SET NOT NULL;
ALTER TABLE fusion ALTER COLUMN strand2_plus SET NOT NULL;
