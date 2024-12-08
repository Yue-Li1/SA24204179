#' @title ADMM Algorithm for Constrained Optimization
#' @description Find the minimum value of f(x)+g(z) under the condition Ax+Bz=C 
#' @param A A numeric matrix.
#' @param B A numeric matrix.
#' @param C A numeric matrix.
#' @param f A equation for X.
#' @param g A equation for X.
#' @param rho A numeric scalar representing the penalty parameter for ADMM (default is 1).
#' @param max_iter An integer representing the maximum number of iterations (default is 100).
#' @param tol A numeric scalar representing the tolerance for convergence (default is 1e-6).
#' @return A list containing:
#'   \item{x}{The optimized \eqn{x} values.}
#'   \item{z}{The optimized \eqn{z} values.}
#'   \item{lambda}{The final Lagrange multipliers.}
#' @examples
#' \dontrun{
#' A <- matrix(c(1, 2, 3, 4), 2, 2)
#' B <- diag(2)
#' C <- c(1, 1)
#' f <- function(x) sum(x^2)
#' g <- function(z) sum(abs(z))
#' result <- admm_algorithm(A, B, C, f, g)
#' print(result$x)
#' }
#' @export
admm_algorithm <- function(A, B, C, f, g, rho = 1, max_iter = 100, tol = 1e-6) {
  # 初始化变量
  # x: 初始值为零向量，大小为 A 的列数
  x <- matrix(0, nrow = ncol(A), ncol = 1)
  # z: 初始值为零向量，大小为 B 的列数
  z <- matrix(0, nrow = ncol(B), ncol = 1)
  # y: 初始值为零向量，大小为 A 的行数
  y <- matrix(0, nrow = nrow(A), ncol = 1)
  # lambda: 初始化拉格朗日乘子为零向量
  lambda <- matrix(0, nrow = nrow(A), ncol = 1)
  
  # 开始ADMM主循环
  for (iter in 1:max_iter) {
    # 更新 x
    # 求解子问题的最优解，利用线性方程组的解析解
    x <- solve(t(A) %*% A + rho * diag(nrow(A)), 
               t(A) %*% (C - B %*% z) + rho * (y - lambda / rho))
    # 更新 z
    # 求解关于 z 的子问题
    z <- solve(t(B) %*% B + rho * diag(nrow(B)), 
               t(B) %*% (C - A %*% x) + rho * (y - lambda / rho))
    # 更新 y
    # 在此问题中，假设 y 的更新直接等于 x
    y <- x
    # 更新拉格朗日乘子 lambda
    # 使用标准的拉格朗日乘子更新公式
    lambda <- lambda + rho * (A %*% x + B %*% z - C)
    # 收敛判断
    # 基于原始残差的平方和判定收敛
    if (sum((A %*% x + B %*% z - C)^2) < tol) {
      cat("Converged at iteration", iter, "\n")
      break
    }
  }
  # 返回优化结果
  return(list(x = x, z = z, lambda = lambda))
}