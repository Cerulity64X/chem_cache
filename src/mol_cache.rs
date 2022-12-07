use std::{collections::HashMap, hash::Hash, error::Error};

use pubchem::{Compound, model::rest::Properties, CompoundProperty};
use serde_json::{value::Serializer, Value, Map, json};

type Prop = CompoundProperty;

const ALL_PROPERTIES: &[CompoundProperty] = &[
    // big property
    CompoundProperty::MolecularFormula,
    CompoundProperty::MolecularWeight,
    CompoundProperty::CanonicalSMILES,
    CompoundProperty::IsomericSMILES,
    CompoundProperty::InChI,
    CompoundProperty::InChIKey,
    CompoundProperty::IUPACName,
    CompoundProperty::Title,
    CompoundProperty::XLogP,
    CompoundProperty::ExactMass,
    CompoundProperty::MonoisotopicMass,
    CompoundProperty::TPSA,
    CompoundProperty::Complexity,
    CompoundProperty::Charge,
    CompoundProperty::HBondDonorCount,
    CompoundProperty::HBondAcceptorCount,
    CompoundProperty::RotatableBondCount,
    CompoundProperty::HeavyAtomCount,
    CompoundProperty::IsotopeAtomCount,
    CompoundProperty::AtomStereoCount,
    CompoundProperty::DefinedAtomStereoCount,
    CompoundProperty::UndefinedAtomStereoCount,
    CompoundProperty::BondStereoCount,
    CompoundProperty::DefinedBondStereoCount,
    CompoundProperty::UndefinedBondStereoCount,
    CompoundProperty::CovalentUnitCount,
    CompoundProperty::Volume3D,
    CompoundProperty::XStericQuadrupole3D,
    CompoundProperty::YStericQuadrupole3D,
    CompoundProperty::ZStericQuadrupole3D,
    CompoundProperty::FeatureCount3D,
    CompoundProperty::FeatureAcceptorCount3D,
    CompoundProperty::FeatureDonorCount3D,
    CompoundProperty::FeatureAnionCount3D,
    CompoundProperty::FeatureCationCount3D,
    CompoundProperty::FeatureRingCount3D,
    CompoundProperty::FeatureHydrophobeCount3D,
    CompoundProperty::ConformerModelRMSD3D,
    CompoundProperty::EffectiveRotorCount3D,
    CompoundProperty::ConformerCount3D,
    CompoundProperty::Fingerprint2D
];

#[derive(Clone, Hash, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct SerCompound {
    pub namespace: String,
    pub identifier: String
}
impl SerCompound {
    // Visibility issues forbid this
    /*pub fn from(cmp: &Compound) -> SerCompound {
        SerCompound { namespace: cmp.namespace.to_string(), identifier: cmp.identifier.to_string() }
    }*/
    pub fn new(id: u32) -> Self {
        Self {
            namespace: String::from("cid"),
            identifier: id.to_string(),
        }
    }
    pub fn with_name(name: &str) -> Self {
        Self {
            namespace: String::from("name"),
            identifier: name.to_string(),
        }
    }
    pub fn with_smiles(smiles: &str) -> Self {
        Self {
            namespace: String::from("smiles"),
            identifier: smiles.to_string(),
        }
    }
    pub fn with_inchi(inchi: &str) -> Self {
        Self {
            namespace: String::from("inchi"),
            identifier: inchi.to_string(),
        }
    }
    pub fn with_inchikey(inchikey: &str) -> Self {
        Self {
            namespace: String::from("inchikey"),
            identifier: inchikey.to_string(),
        }
    }
    pub fn to_compound(&self) -> Result<Option<Compound>, Box<dyn Error>> {
        match &self.namespace[..] {
            "cid" => Ok(Some(Compound::new(self.identifier.parse::<u32>()?))),
            "name" => Ok(Some(Compound::with_name(&self.identifier))),
            "smiles" => Ok(Some(Compound::with_smiles(&self.identifier))),
            "inchi" => Ok(Some(Compound::with_inchi(&self.identifier))),
            "inchikey" => Ok(Some(Compound::with_inchikey(&self.identifier))),
            _ => Ok(None)
        }
    }
}

#[derive(Debug)]
pub struct CompoundCache {
    cache: HashMap<SerCompound, Properties>
}
impl CompoundCache {
    pub fn new() -> CompoundCache {
        CompoundCache { cache: HashMap::new() }
    }
    /// Use overwrite for overwriting, this will not insert if value exists. If the compound namespaces are not the same, then the compound properties will be duplicated.
    pub fn store(&mut self, cmp: SerCompound) -> Result<(), Box<dyn Error>> {
        if !self.cache.contains_key(&cmp) {
            let props = cmp.to_compound()?.ok_or(String::new())?.properties(ALL_PROPERTIES)?;
            self.cache.insert(cmp, props);
        }
        Ok(())
    }
    /// Overwrites properties.
    pub fn overwrite(&mut self, cmp: SerCompound) -> Result<(), Box<dyn Error>>{
        let props = cmp.to_compound()?.ok_or(String::new())?.properties(ALL_PROPERTIES)?;
        self.cache.insert(cmp, props);
        Ok(())
    }
    /// If the compound does not exist, the properties are added and returned.
    pub fn get(&mut self, cmp: SerCompound) -> Result<(bool, &Properties), Box<dyn Error>> {
        let props = cmp.to_compound()?.ok_or(String::new())?.properties(ALL_PROPERTIES)?;
        let haskey = self.cache.contains_key(&cmp);
        if !haskey {
            self.cache.insert(cmp.clone(), props);
        }
        Ok((haskey, &self.cache[&cmp]))
    }
    /// If the compound does not exist, None is returned. Does not make a PubChem request.
    pub fn get_noreq(&self, cmp: SerCompound) -> Result<Option<&Properties>, pubchem::error::Error> {
        if self.cache.contains_key(&cmp) {
            Ok(Some(&self.cache[&cmp]))
        } else {
            Ok(None)
        }
    }

    pub fn insert(&mut self, key: SerCompound, val: Properties) {
        self.cache.insert(key, val);
    }

    pub fn serialize(&self) -> Result<Value, String> {
        let mut arr: Vec<Value> = Vec::new();
        for (cmp, prop) in &self.cache {
            let mut ser_obj = Map::new();
            ser_obj.insert("namespace".to_owned(), Value::String(cmp.namespace.clone()));
            ser_obj.insert("identifier".to_owned(), Value::String(cmp.identifier.clone()));
            let mut properties = Map::new();
            {
                // big property 2 electrig boogaloo
                properties.insert("atom_stereo_count".to_owned(), prop.atom_stereo_count.unwrap().into());
                properties.insert("bond_stereo_count".to_owned(), prop.bond_stereo_count.unwrap().into());
                properties.insert("canonical_smiles".to_owned(), Value::String(prop.canonical_smiles.as_ref().unwrap().clone()));
                properties.insert("charge".to_owned(), prop.charge.unwrap().into());
                properties.insert("cid".to_owned(), prop.cid.into());
                properties.insert("complexity".to_owned(), prop.complexity.unwrap().into());
                properties.insert("conformer_count_3d".to_owned(), prop.conformer_count_3d.unwrap().into());
                properties.insert("conformer_model_rmsd_3d".to_owned(), prop.conformer_model_rmsd_3d.unwrap().into());
                properties.insert("covalent_unit_count".to_owned(), prop.covalent_unit_count.unwrap().into());
                properties.insert("defined_atom_stereo_count".to_owned(), prop.defined_atom_stereo_count.unwrap().into());
                properties.insert("defined_bond_stereo_count".to_owned(), prop.defined_bond_stereo_count.unwrap().into());
                properties.insert("effective_rotor_count_3d".to_owned(), prop.effective_rotor_count_3d.unwrap().into());
                properties.insert("exact_mass".to_owned(), Value::String(prop.exact_mass.as_ref().unwrap().clone()));
                properties.insert("feature_acceptor_count_3d".to_owned(), prop.feature_acceptor_count_3d.unwrap().into());
                properties.insert("feature_anion_count_3d".to_owned(), prop.feature_anion_count_3d.unwrap().into());
                properties.insert("feature_cation_count_3d".to_owned(), prop.feature_cation_count_3d.unwrap().into());
                properties.insert("feature_count_3d".to_owned(), prop.feature_count_3d.unwrap().into());
                properties.insert("feature_donor_count_3d".to_owned(), prop.feature_donor_count_3d.unwrap().into());
                properties.insert("feature_hydrophobe_count_3d".to_owned(), prop.feature_hydrophobe_count_3d.unwrap().into());
                properties.insert("feature_ring_count_3d".to_owned(), prop.feature_ring_count_3d.unwrap().into());
                properties.insert("fingerprint_2d".to_owned(), Value::String(prop.fingerprint_2d.as_ref().unwrap().clone()));
                properties.insert("hbond_acceptor_count".to_owned(), prop.hbond_acceptor_count.unwrap().into());
                properties.insert("hbond_donor_count".to_owned(), prop.hbond_donor_count.unwrap().into());
                properties.insert("heavy_atom_count".to_owned(), prop.heavy_atom_count.unwrap().into());
                properties.insert("inchi".to_owned(), Value::String(prop.inchi.as_ref().unwrap().clone()));
                properties.insert("inchi_key".to_owned(), Value::String(prop.inchi_key.as_ref().unwrap().clone()));
                properties.insert("isomeric_smiles".to_owned(), Value::String(prop.isomeric_smiles.as_ref().unwrap().clone()));
                properties.insert("isotope_atom_count".to_owned(), prop.isotope_atom_count.unwrap().into());
                properties.insert("iupac_name".to_owned(), valify_string(&prop.iupac_name));
                properties.insert("molecular_formula".to_owned(), Value::String(prop.molecular_formula.as_ref().unwrap().clone()));
                properties.insert("molecular_weight".to_owned(), Value::String(prop.molecular_weight.as_ref().unwrap().clone()));
                properties.insert("monoisotopic_mass".to_owned(), Value::String(prop.monoisotopic_mass.as_ref().unwrap().clone()));
                properties.insert("rotatable_bond_count".to_owned(), prop.rotatable_bond_count.unwrap().into());
                properties.insert("title".to_owned(), Value::String(prop.title.as_ref().unwrap().clone()));
                properties.insert("tpsa".to_owned(), prop.tpsa.unwrap().into());
                properties.insert("undefined_atom_stereo_count".to_owned(), prop.undefined_atom_stereo_count.unwrap().into());
                properties.insert("undefined_bond_stereo_count".to_owned(), prop.undefined_bond_stereo_count.unwrap().into());
                properties.insert("volume_3d".to_owned(), prop.volume_3d.unwrap().into());
                properties.insert("x_steric_quadrupole_3d".to_owned(), prop.x_steric_quadrupole_3d.unwrap().into());
                properties.insert("xlogp".to_owned(), prop.xlogp.unwrap().into());
                properties.insert("y_steric_quadrupole_3d".to_owned(), prop.y_steric_quadrupole_3d.unwrap().into());
                properties.insert("z_steric_quadrupole_3d".to_owned(), prop.z_steric_quadrupole_3d.unwrap().into());
            }
            ser_obj.insert("properties".to_owned(), Value::Object(properties));
            arr.push(Value::Object(ser_obj));
        }
        let mut map = Map::new();
        map.insert("cache".to_owned(), Value::Array(arr));
        Ok(Value::Object(map))
    }

    pub fn deserialize(st: String) -> Result<CompoundCache, String> {
        let mut output_cache = CompoundCache::new();
        match serde_json::from_str::<Value>(&st[..]) {
            Ok(root) => {
                let cache = root
                    .as_object().ok_or("The root JSON was not an object!")?
                    .get("cache").ok_or("`cache` could not be found! Make sure it's an array in the root object!")?
                    .as_array().ok_or("`cache` was not an array!")?;
                for i in cache {
                    match i {
                        Value::Object(entry) => {
                            let obj = &i["properties"];
                            let properties = Properties {
                                // big property 3: deser
                                atom_stereo_count: Some(obj["atom_stereo_count"].as_i64().unwrap() as i32),
                                bond_stereo_count: Some(obj["bond_stereo_count"].as_i64().unwrap() as i32),
                                canonical_smiles: Some(obj["canonical_smiles"].as_str().clone().unwrap().to_owned()),
                                charge: Some(obj["charge"].as_i64().unwrap() as i32),
                                cid: obj["cid"].as_i64().unwrap() as i32,
                                complexity: Some(obj["complexity"].as_i64().unwrap() as i32),
                                conformer_count_3d: Some(obj["conformer_count_3d"].as_i64().unwrap() as i32),
                                conformer_model_rmsd_3d: Some(obj["conformer_model_rmsd_3d"].as_f64().unwrap()),
                                covalent_unit_count: Some(obj["covalent_unit_count"].as_i64().unwrap() as i32),
                                defined_atom_stereo_count: Some(obj["defined_atom_stereo_count"].as_i64().unwrap() as i32),
                                defined_bond_stereo_count: Some(obj["defined_bond_stereo_count"].as_i64().unwrap() as i32),
                                effective_rotor_count_3d: Some(obj["effective_rotor_count_3d"].as_f64().unwrap()),
                                exact_mass: Some(obj["exact_mass"].as_str().unwrap().to_owned()),
                                feature_acceptor_count_3d: Some(obj["feature_acceptor_count_3d"].as_i64().unwrap() as i32),
                                feature_anion_count_3d: Some(obj["feature_anion_count_3d"].as_i64().unwrap() as i32),
                                feature_cation_count_3d: Some(obj["feature_cation_count_3d"].as_i64().unwrap() as i32),
                                feature_count_3d: Some(obj["feature_count_3d"].as_i64().unwrap() as i32),
                                feature_donor_count_3d: Some(obj["feature_donor_count_3d"].as_i64().unwrap() as i32),
                                feature_hydrophobe_count_3d: Some(obj["feature_hydrophobe_count_3d"].as_i64().unwrap() as i32),
                                feature_ring_count_3d: Some(obj["feature_ring_count_3d"].as_i64().unwrap() as i32),
                                fingerprint_2d: Some(obj["fingerprint_2d"].as_str().unwrap().to_owned()),
                                hbond_acceptor_count: Some(obj["hbond_acceptor_count"].as_i64().unwrap() as i32),
                                hbond_donor_count: Some(obj["hbond_donor_count"].as_i64().unwrap() as i32),
                                heavy_atom_count: Some(obj["heavy_atom_count"].as_i64().unwrap() as i32),
                                inchi: Some(obj["inchi"].as_str().unwrap().to_owned()),
                                inchi_key: Some(obj["inchi_key"].as_str().unwrap().to_owned()),
                                isomeric_smiles: Some(obj["isomeric_smiles"].as_str().unwrap().to_owned()),
                                isotope_atom_count: Some(obj["isotope_atom_count"].as_i64().unwrap() as i32),
                                iupac_name: obj["iupac_name"].as_str().map(|s|s.to_owned()),
                                molecular_formula: Some(obj["molecular_formula"].as_str().unwrap().to_owned()),
                                molecular_weight: Some(obj["molecular_weight"].as_str().unwrap().to_owned()),
                                monoisotopic_mass: Some(obj["monoisotopic_mass"].as_str().unwrap().to_owned()),
                                rotatable_bond_count: Some(obj["rotatable_bond_count"].as_i64().unwrap() as i32),
                                title: Some(obj["title"].as_str().unwrap().to_owned()),
                                tpsa: Some(obj["tpsa"].as_f64().unwrap()),
                                undefined_atom_stereo_count: Some(obj["undefined_atom_stereo_count"].as_i64().unwrap() as i32),
                                undefined_bond_stereo_count: Some(obj["undefined_bond_stereo_count"].as_i64().unwrap() as i32),
                                volume_3d: Some(obj["volume_3d"].as_f64().unwrap()),
                                x_steric_quadrupole_3d: Some(obj["x_steric_quadrupole_3d"].as_f64().unwrap()),
                                xlogp: Some(obj["xlogp"].as_f64().unwrap()),
                                y_steric_quadrupole_3d: Some(obj["y_steric_quadrupole_3d"].as_f64().unwrap()),
                                z_steric_quadrupole_3d: Some(obj["z_steric_quadrupole_3d"].as_f64().unwrap())
                            };
                            let key = SerCompound {
                                namespace: entry["namespace"].as_str().unwrap().to_owned(),
                                identifier: entry["identifier"].as_str().unwrap().to_owned()
                            };
                            output_cache.insert(key, properties);
                        }
                        _ => Err("Value was not an object!")?
                    }
                }
            }
            Err(e) => Err(format!("Could not parse JSON! ({e})"))?
        }
        Ok(output_cache)
    }
}

pub fn valify_string(string: &Option<String>) -> Value {
    match string {
        Some(st) => Value::String(st.clone()),
        None => Value::Null
    }
}
