use bgzip::BGZFReader;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use serde_json;
use std::char;
use std::fmt;
use std::fs::File;
use std::io::BufRead;
use std::path::Path;
use std::rc::Rc;

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct GFFRecord {
    pub seqid: String,
    pub source: String,
    #[serde(rename = "type")]
    pub type_name: String,
    pub bgn: u32,
    pub end: u32,
    pub score: Option<f32>,
    pub strand: char,
    pub phase: Option<u8>,
    pub attributes: FxHashMap<String, String>,
}

impl GFFRecord {
    pub fn from_line(line: &String) -> GFFRecord {
        let fields = line
            .trim_end()
            .split("\t")
            .into_iter()
            .map(|s| s.to_string())
            .collect::<Vec<String>>();
        GFFRecord::from_fields(&fields)
    }

    pub fn from_fields(fields: &Vec<String>) -> GFFRecord {
        let seqid = fields[0].clone();
        let source = fields[1].clone();
        let type_ = fields[2].clone();
        let bgn = fields[3]
            .parse::<u32>()
            .expect(format!("can't parse the start coordindate").as_str());
        let end = fields[4]
            .parse::<u32>()
            .expect(format!("can't parse the end coordindate").as_str());
        let score = match fields[5].as_str() {
            "." => None,
            s => Some(
                s.parse::<f32>()
                    .expect(format!("can't parse score").as_str()),
            ),
        };

        let strand = fields[6][0..1].chars().nth(0).unwrap();
        let phase =
            match fields[7].as_str() {
                "." => None,
                s => Some(s.parse::<u8>().expect(
                    format!("fail to parse the phase field {}", fields[6].as_str(),).as_str(),
                )),
            };
        let attributes = fields[8]
            .split(";")
            .into_iter()
            .map(|s| {
                let kv = s
                    .split("=")
                    .into_iter()
                    .map(|s| s.to_string())
                    .collect::<Vec<String>>();
                if kv.len() != 2 {
                    panic!("error parsing attributes")
                };
                (kv[0].clone(), kv[1].clone())
            })
            .collect::<FxHashMap<String, String>>();

        Self {
            seqid,
            source,
            type_name: type_,
            bgn,
            end,
            score,
            strand,
            phase,
            attributes,
        }
    }
}

impl fmt::Display for GFFRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut out = vec![];
        out.push(format!(
            "{}\t{}\t{}\t{}\t{}",
            self.seqid, self.source, self.type_name, self.bgn, self.end
        ));

        out.push(if self.score.is_none() {
            ".".to_string()
        } else {
            format!("{}", self.score.unwrap())
        });

        out.push(format!("{}", self.strand));

        out.push(if self.phase.is_none() {
            ".".to_string()
        } else {
            format!("{}", self.phase.unwrap())
        });

        out.push(
            self.attributes
                .iter()
                .map(|(k, v)| format!("{}={}", k, v))
                .collect::<Vec<String>>()
                .join(";"),
        );

        write!(f, "{}", out.join("\t"))
    }
}

type IdToGffRec = FxHashMap<String, Rc<GFFRecord>>;
type IdToChildren = FxHashMap<String, Vec<Rc<GFFRecord>>>;
type NameToGffRec = FxHashMap<String, Rc<GFFRecord>>;

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct GFFDB {
    pub header: Vec<String>,
    pub records: Vec<Rc<GFFRecord>>,
    pub id_to_rec: IdToGffRec,
    pub name_to_rec: NameToGffRec,
    pub children: IdToChildren,
}

impl GFFDB {
    pub fn from_bgzip_file(filepath: &Path) -> std::io::Result<GFFDB> {
        let file = BGZFReader::new(File::open(filepath)?);
        let mut header = Vec::<String>::new();
        let mut records = Vec::<Rc<GFFRecord>>::new();
        let mut id_to_rec = IdToGffRec::default();
        let mut name_to_rec = NameToGffRec::default();
        let mut children = IdToChildren::default();
        file.lines().into_iter().for_each(|line| {
            let line = line.unwrap();
            if &line[0..1] != "#" {
                let rec = Rc::new(GFFRecord::from_line(&line));
                records.push(rec.clone());

                if rec.attributes.contains_key("ID") {
                    let id = rec.attributes.get("ID").unwrap();
                    id_to_rec.insert(id.clone(), rec.clone());
                }
                if rec.attributes.contains_key("Name") {
                    let name = rec.attributes.get("Name").unwrap();
                    name_to_rec.insert(name.clone(), rec.clone());
                }
                if rec.attributes.contains_key("Parent") {
                    let parent_id = rec.attributes.get("Parent").unwrap();
                    children
                        .entry(parent_id.clone())
                        .or_insert_with(|| vec![])
                        .push(rec.clone());
                }
            } else {
                header.push(line);
            }
        });
        Ok(GFFDB {
            header,
            records,
            id_to_rec,
            name_to_rec,
            children,
        })
    }

    pub fn from_list_of_fields(list_of_fields: &Vec<Vec<String>>) -> GFFDB {
        let header = Vec::<String>::new();
        let mut records = Vec::<Rc<GFFRecord>>::new();
        let mut id_to_rec = IdToGffRec::default();
        let mut name_to_rec = NameToGffRec::default();
        let mut children = IdToChildren::default();

        list_of_fields.into_iter().for_each(|fields| {
            let rec = Rc::new(GFFRecord::from_fields(&fields));
            records.push(rec.clone());

            if rec.attributes.contains_key("ID") {
                let id = rec.attributes.get("ID").unwrap();
                id_to_rec.insert(id.clone(), rec.clone());
            }
            if rec.attributes.contains_key("Name") {
                let name = rec.attributes.get("Name").unwrap();
                name_to_rec.insert(name.clone(), rec.clone());
            }
            if rec.attributes.contains_key("Parent") {
                let parent_id = rec.attributes.get("Parent").unwrap();
                children
                    .entry(parent_id.clone())
                    .or_insert_with(|| vec![])
                    .push(rec.clone());
            }
        });

        GFFDB {
            header,
            records,
            id_to_rec,
            name_to_rec,
            children,
        }
    }

    pub fn get_all_offspring(
        &self,
        id_or_name: &String,
        recusive: bool,
    ) -> Option<Vec<Rc<GFFRecord>>> {
        let mut all_offspring = Vec::<Rc<GFFRecord>>::new();

        let id = if self.id_to_rec.contains_key(id_or_name) {
            Some(id_or_name)
        } else {
            if self.name_to_rec.contains_key(id_or_name) {
                let r = self.name_to_rec.get(id_or_name).unwrap();
                r.attributes.get("ID")
            } else {
                None
            }
        };

        if id.is_none() {
            return None;
        }

        let id = id.unwrap();
        match self.children.get(id) {
            Some(children) => {
                children.iter().for_each(|r| {
                    if recusive && r.attributes.contains_key("ID") {
                        let id = r.attributes.get("ID").unwrap();
                        if let Some(more_offsprings) = self.get_all_offspring(id, recusive) {
                            more_offsprings.iter().for_each(|r| {
                                all_offspring.push(r.clone());
                            });
                        }
                    }
                    all_offspring.push(r.clone());
                });
                Some(all_offspring)
            }
            None => None,
        }
    }

    pub fn dump_json(&self) {
        println!("{}", serde_json::to_string(&self).unwrap());
    }

    pub fn load_json(s: &String) -> serde_json::Result<GFFDB> {
        let gffdb: GFFDB = serde_json::from_str(s)?;
        Ok(gffdb)
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct QueryOut {
    parent: Rc<GFFRecord>,
    offspring: Vec<Rc<GFFRecord>>,
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_gff_to_db() {
        let res = super::GFFDB::from_bgzip_file(&Path::new("./test/test_data/test.gff3.gz"));
        let gdb = res.unwrap();
        println!("{}", gdb.header.join("\n"));
        let r = gdb.name_to_rec.get(&"FLG".to_string()).unwrap();
        let parent = r.clone();
        let mut offspring = Vec::<Rc<GFFRecord>>::new();
        println!("{}", r);
        gdb.get_all_offspring(&"FLG".to_string(), true)
            .unwrap()
            .into_iter()
            .for_each(|r| {
                println!("{}", r);
                offspring.push(r.clone());
            });
        let qr = QueryOut { parent, offspring };
        println!("{}", serde_json::to_string(&qr).unwrap());
    }
}
