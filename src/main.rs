use std::{error::Error, env::args, fs};

use pubchem::Compound;

use crate::mol_cache::{CompoundCache, SerCompound};

mod mol_cache;

fn main() -> Result<(), Box<dyn Error>> {
    // load compounds
    let mut cache = CompoundCache::deserialize(
        fs::read_to_string("compounds.json").unwrap_or(String::from("{}"))
    ).unwrap_or(CompoundCache::new());

    cache.get(SerCompound::with_smiles("O"))?;

    // write compounds
    fs::write("compounds.json", cache.serialize()?.to_string())?;
    println!("Wrote to compounds.json.");
    Ok(())
}
