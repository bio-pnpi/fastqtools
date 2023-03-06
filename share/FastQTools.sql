DROP TABLE IF EXISTS header; 
CREATE TABLE header (
  id UUID PRIMARY KEY,
  runid TEXT,
  sampleid TEXT,
  read INT,
  ch INT,
  start_time TEXT,
  model_version_id TEXT,
  barcode TEXT);

DROP TABLE IF EXISTS sequence;
CREATE TABLE sequence (
	id UUID PRIMARY KEY,
	sequence TEXT,
	quality TEXT,
	len INT
);

DROP TABLE IF EXISTS dprimers;
CREATE TABLE dprimers (
	id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
	primer TEXT NOT NULL UNIQUE,
	len INT
);

DROP TABLE IF EXISTS fprimers;
CREATE TABLE fprimers (
	id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
	primer TEXT NOT NULL,
	dprimer UUID NOT NULL,
	UNIQUE(primer, dprimer)
);

DROP TABLE IF EXISTS rprimers;
CREATE TABLE rprimers (
        id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
        primer TEXT NOT NULL,
        dprimer UUID NOT NULL,
        UNIQUE(primer, dprimer)
);
