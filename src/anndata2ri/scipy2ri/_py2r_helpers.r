sparse_matrix <- function(x, conv_data, dims, ...) {
    Matrix::sparseMatrix(
        ...,
        x=conv_data(x),
        dims=as.integer(dims),
        index1=FALSE
    )
}

from_csc <- function(i, p, x, dims, conv_data) {
    sparse_matrix(
        i=as.integer(i),
        p=as.integer(p),
        x=x,
        conv_data=conv_data,
        dims=dims,
        repr="C"
    )
}

from_csr <- function(j, p, x, dims, conv_data) {
    sparse_matrix(
        j=as.integer(j),
        p=as.integer(p),
        x=x,
        conv_data=conv_data,
        dims=dims,
        repr="R"
    )
}

from_coo <- function(i, j, x, dims, conv_data) {
    sparse_matrix(
        i=as.integer(i),
        j=as.integer(j),
        x=x,
        conv_data=conv_data,
        dims=dims,
        repr="T"
    )
}

from_dia <- function(n, x, conv_data) {
    Matrix::Diagonal(n=as.integer(n), x=conv_data(x))
}
