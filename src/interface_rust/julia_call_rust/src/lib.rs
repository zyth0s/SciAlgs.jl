#[no_mangle]
pub extern fn double_input(input: i32) -> i32 {
    println!("Hello from Rust");
    input * 2
}

#[cfg(test)]
mod tests {
    use crate::double_input;
    #[test]
    fn it_works() {
        //assert_eq!(2 + 2, 4);
        assert_eq!(double_input(2), 4);
    }
}
