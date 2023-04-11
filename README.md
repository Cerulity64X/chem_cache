# Chemical Storage
Crate containing a serializable storage for the PubChem database.\
The `pubchem` crate frequently makes requests, this crate provides a means of storing those requests to avoid repeats.\
The storage is serializable/deserializable to JSON, meaning you can save/load compounds easily.
# Usage
`CompoundCache` is a struct that contains information for multiple chemicals. It has functions for loading, storing, getting, and saving elements. `SerCompound` is a struct that defines a queryable element. This is also the key type for `CompoundCache`.
# Plans
- Serde to more compact formats
- Time serialization and expiration
