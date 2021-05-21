extern crate ndarray;
use ndarray::prelude::*;

#[no_mangle]
pub extern fn double_input(input: i32) -> i32 {
    println!("Hello from Rust!");
    input * 2
}

// Take a pointer to the first element of a C array
// with len elements of type f64, assuming some properties (unsafe)
#[no_mangle]
pub extern fn modify_vectorf64(v_ptr: *mut f64, len: i64) {
    let mut v = unsafe { ArrayViewMut::from_shape_ptr(len as usize, v_ptr) };
    v.fill(3.3);
    // Indexing starts at 0; assuming len>1 !!
    v[[1]] = 4.5;
    println!("Your vector has been modified: v = {}", v);
}

// Take a pointer to the first element of a C array
// with nrow*ncol elements of type f64, assuming some properties (unsafe)
#[no_mangle]
pub extern fn modify_matrixf64(m_ptr: *mut f64, nrow: i64, ncol: i64) {
    // - Default   :    row-major layout (like C)
    // - shape.f() : column-major layout (like Fortran)
    let mut m = unsafe { ArrayViewMut::from_shape_ptr((nrow as usize, ncol as usize).f(), m_ptr) };
    m.fill(0.0);
    // Indexing starts at 0; assuming nrow>1 & ncol>1 !!
    m[[1,1]] = 4.5;
    println!("Your matrix has been modified: m =");
    println!("{}", m);
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_double_input() {
        assert_eq!(double_input(2), 4);
    }

    #[test]
    fn test_modify_vectorf64() {
        let mut v = array![0.0,0.0,0.0];
        modify_vectorf64(v.as_mut_ptr(), 3);
        assert_eq!(v, array![3.3, 4.5, 3.3]);
    }

    #[test]
    fn test_modify_matrixf64() {
        let mut m = array![[0.0, 0.0], [0.0, 0.0]];
        modify_matrixf64(m.as_mut_ptr(), 2, 2);
        let m_ref = array![[0.0, 0.0], [0.0, 4.5]];
        assert_eq!(m, m_ref);
    }
}
