/// Include the generated bindings into a separate module.
#[allow(non_upper_case_globals)]
#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
#[allow(unused)]
pub mod wfa {
    include!(concat!(env!("OUT_DIR"), "/bindings_wfa.rs"));
}
