use super::{
    arithmetic::{add2, sub2, div_rem},
    Field, Zero, One,
};

use core::{
    debug_assert, unimplemented,
    clone::Clone,
    cmp::{PartialEq, PartialOrd, Eq, Ord, Ordering},
    marker::{Copy, PhantomData},
};

pub trait FpConfig<const N: usize> {
    const MODULUS: [u64; N];
    const ONE: [u64; N];  // Same as Montgomery Radix
    const ZERO: [u64; N] = [0; N];
}

pub trait PrimeField<P: FpConfig<N>, const N: usize>: Field<N> {
    fn reduce(&mut self);
    // fn to_uint(&self);
}

pub struct Fp<P: FpConfig<N>, const N: usize> {
    pub repr: [u64; N],
    pub phantom: PhantomData<P>
}

impl<P: FpConfig<N>, const N: usize> Clone for Fp<P, N> {
    fn clone(&self) -> Fp<P, N> {
        Self {repr: self.repr.clone(), phantom: PhantomData}
    }
}

impl<P: FpConfig<N>, const N: usize> Copy for Fp<P, N> {}

impl<P: FpConfig<N>, const N: usize> PartialEq for Fp<P, N> {
    fn eq(&self, other: &Self) -> bool {
        self.repr == other.repr
    }

    fn ne(&self, other: &Self) -> bool {
        self.repr != other.repr
    }
}

impl<P: FpConfig<N>, const N: usize> Eq for Fp<P, N> {}

impl<P: FpConfig<N>, const N: usize> PartialOrd for Fp<P, N> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.repr.partial_cmp(&other.repr)
    }
}

impl<P: FpConfig<N>, const N: usize> Ord for Fp<P, N> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.repr.cmp(&other.repr)
    }
}

impl<P: FpConfig<N>, const N: usize> PrimeField<P, N> for Fp<P, N> {
    fn reduce(&mut self) {
        div_rem(&mut self.repr, &P::MODULUS, 0);
    }
}

impl<P: FpConfig<N>, const N: usize> Default for Fp<P, N> {
    fn default() -> Self {
        Self::zero()
    }
}

impl<P: FpConfig<N>, const N: usize> Zero for Fp<P, N> {
    fn zero() -> Self {
        Self {repr: [0u64; N], phantom: PhantomData}
    }

    fn is_zero(&self) -> bool {
        self.repr == Self::zero().repr
    }
}

impl<P: FpConfig<N>, const N: usize> One for Fp<P, N> {
    fn one() -> Self {
        let mut data = [0u64; N];
        data[0] = 1;
        Self {repr: data, phantom: PhantomData}
    }

    fn is_one(&self) -> bool {
        self.repr == Self::one().repr
    }
}

impl<P: FpConfig<N>, const N: usize> AsRef<[u64; N]> for Fp<P, N> {
    fn as_ref(&self) -> &[u64; N] {
        &self.repr
    }
}

impl<P: FpConfig<N>, const N: usize> AsMut<[u64; N]> for Fp<P, N> {
    fn as_mut(&mut self) -> &mut [u64; N] {
        &mut self.repr
    }
}

impl<P: FpConfig<N>, const N: usize> From<[u64; N]> for Fp<P, N> {
    fn from(value: [u64; N]) -> Self {
        let mut val = value.clone();
        let _ = div_rem(&mut val, &P::MODULUS, 0);
        Self {repr: val, phantom: PhantomData}
    }
}

impl<P: FpConfig<N>, const N: usize> Field<N> for Fp<P, N> {
    fn add_assign(&mut self, other: &Self) {
        let carry = add2(&mut self.repr, &other.repr);
        if carry > 0 || self.repr > P::MODULUS {
            let borrow = sub2(&mut self.repr, &P::MODULUS);
            // debug_assert!(u64::from(borrow) == carry);
        }

    }

    fn mul_assign(&mut self, other: &Self) {
        unimplemented!();

        // Least-significant zeros have no effect on the output.
        // if let Some(&0) = other.0[0] {
        // if let Some(nz) = other.0.iter().position(|&d| d != 0) {
        // b = &b[nz..];
        // acc = &mut acc[nz..];
        // } else {
        // return;
        // }
        // }
        // if let Some(&0) = c.first() {
        // if let Some(nz) = c.iter().position(|&d| d != 0) {
        // c = &c[nz..];
        // acc = &mut acc[nz..];
        // } else {
        // return;
        // }
        // }

        // let acc = acc;
        // let (x, y) = if b.len() < c.len() { (b, c) } else { (c, b) };

        // Long multiplication:
        // for (i, xi) in x.iter().enumerate() {
        //     mac_digit(&mut acc[i..], y, *xi);
        // }
    }

    //    if c == 0 {
    //        return;
    //    }
    //
    //    let mut carry = 0;
    //    let (a_lo, a_hi) = acc.split_at_mut(b.len());
    //
    //    for (a, &b) in a_lo.iter_mut().zip(b) {
    //        *a = mac_with_carry(*a, b, c, &mut carry);
    //    }
    //
    //    let (carry_hi, carry_lo) = big_digit::from_doublebigdigit(carry);
    //
    //    let final_carry = if carry_hi == 0 {
    //        __add2(a_hi, &[carry_lo])
    //    } else {
    //        __add2(a_hi, &[carry_hi, carry_lo])
    //    };
    //    assert_eq!(final_carry, 0, "carry overflow during multiplication!");

    fn pow_assign(&mut self, exp: u32) {
        if exp == 0 {
            *self = Self::one();
            return ();
        }
        if *self == Self::zero() || *self == Self::one() {
            return ();
        }
        let mut i = exp;
        let mut acc = self.clone();
        *self = Self::one();
        while i != 0 { 
            if i % 2 == 1 {
                self.mul_assign(&acc);
            }
            acc.mul_assign(&acc.clone());
            i >>= 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::{
        clone::Clone,
        marker::Copy,
    };

    pub struct P;
    impl FpConfig<2> for P {
        const MODULUS: [u64; 2] = [5, 1];
        const ONE: [u64; 2] = [1, 1];
    }

    type Fp128 = Fp::<P, 2>;

    #[test]
    fn from_no_division() {
        let input: [u64; 2] = [2, 1];
        let output = Fp128::from(input);
        assert_eq!(input, output.repr);
    }

    #[test]
    fn from_with_division() {
        let input: [u64; 2] = [20, 3];
        let output = Fp128::from(input);
        assert_eq!([5, 0], output.repr);
    }

    #[test]
    fn add_assign_no_carry() {
        let mut input = Fp128::from([2, 1]);
        let other = Fp128::from([2, 0]);
        input.add_assign(&other);
        assert_eq!(*input.as_ref(), [4, 1]);
    }

    #[test]
    fn add_assign_with_carry() {
        let mut input = Fp128::from([(2u128.pow(64) - 1) as u64, 1]);
        let other = Fp128::from([12, 0]);
        input.add_assign(&other);
        assert_eq!(*input.as_ref(), [1, 0]);
    }

    #[test]
    fn pow_assign_zero_exp() {
        let mut input = Fp128::from([4, 3]);
        input.pow_assign(0);
        assert_eq!(*input.as_ref(), [1, 0]);
    }
}
