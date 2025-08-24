MODULE mecanica
  IMPLICIT NONE
CONTAINS

FUNCTION energia_potencial (G, ms, qs, soft) RESULT(V)
  REAL(8), INTENT(IN) :: G, ms(:)
  REAL(8), INTENT(IN) :: qs(:,:)
  REAL(8), INTENT(IN), OPTIONAL :: soft
  REAL(8) :: eps, V, dist2
  INTEGER :: a, b
  ! Amortecimento
  eps = 0.0
  IF (PRESENT(soft)) eps = soft

  V = 0.0
  DO a = 2, SIZE(qs,1)
      DO b = 1, a - 1
          dist2 = DOT_PRODUCT(qs(a,:) - qs(b,:), qs(a,:) - qs(b,:))
          V = V - ms(a) * ms(b) / SQRT(dist2 + eps*eps)
      END DO
  END DO
  V = G * V
END FUNCTION energia_potencial

FUNCTION energia_cinetica (m, P) RESULT(ec)
  REAL(8), INTENT(IN) :: m(:), P(:,:)
  REAL(8) :: ec
  INTEGER  :: i
  ec=0.0
  DO i=1, SIZE(m)
    ec = ec + DOT_PRODUCT(P(i,:), P(i,:))/m(i)
  END DO
  ec = 0.5 * ec
END FUNCTION energia_cinetica

FUNCTION energia_total (G, m, R, P, soft) RESULT(e_tot)
  REAL(8), INTENT(IN) :: m(:), R(:,:), P(:,:), G
  REAL(8), INTENT(IN), OPTIONAL :: soft
  REAL(8) :: e_tot, eps, dist2
  INTEGER :: a, b
  eps = 0.0
  IF (PRESENT(soft)) eps = soft

  e_tot = 0.0
  DO a = 1, SIZE(m)
    e_tot = e_tot + 0.5 * DOT_PRODUCT(P(a,:),P(a,:)) / m(a)
    DO b = 1, a - 1
      dist2 = DOT_PRODUCT(R(a,:) - R(b,:), R(a,:) - R(b,:))
      e_tot = e_tot - G * m(a) * m(b) / SQRT(dist2 + eps*eps)
    END DO
  END DO
END FUNCTION energia_total

FUNCTION energia_total_separada (G, m, R, P, soft) RESULT(e_tots)
  REAL(8), INTENT(IN) :: m(:), R(:,:), P(:,:), G
  REAL(8), OPTIONAL, INTENT(IN) :: soft
  REAL(8) :: e_tots(SIZE(m))
  REAL(8) :: eps2, pot, dist(3)
  INTEGER :: a, b

  eps2 = 0.0
  IF (PRESENT(soft)) eps2 = soft*soft

  e_tots = 0.0

  DO a = 1, SIZE(m)
    ! Contribuicao cinetica
    e_tots(a) = e_tots(a) + DOT_PRODUCT(P(a,:), P(a,:)) / (2.0 * m(a))

    ! Contribuicao potencial
    DO b = 1, a - 1
      dist = R(b,:) - R(a,:)
      pot = - G * m(a) * m(b) / SQRT(DOT_PRODUCT(dist,dist) + eps2)
      e_tots(a) = e_tots(a) + pot
      e_tots(b) = e_tots(b) + pot
    END DO
  END DO
END FUNCTION energia_total_separada

FUNCTION momento_inercia_separado (m, R) RESULT(is)
  REAL(8), INTENT(IN) :: m(:), R(:,:)
  REAL(8) :: is(SIZE(m)), is3d(SIZE(m),3), termo(3)
  INTEGER :: a, b

  is = 0.0
  is3d = 0.0

  DO a = 2, SIZE(m)
    DO b = 1, a - 1
      termo = m(b) * (R(a,:) - R(b,:))
      is3d(a,:) = is3d(a,:) + termo
      is3d(b,:) = is3d(b,:) + termo
    END DO
  END DO

  DO a = 1, SIZE(m)
    is(a) = m(a) * DOT_PRODUCT(is3d(a,:),is3d(a,:))
  END DO

  is = is / (SUM(m)**2)
END FUNCTION momento_inercia_separado

FUNCTION virial_potencial_amortecido_separado (G, m, R, eps) RESULT(f_prod_q)
  REAL(8), INTENT(IN) :: G, m(:), R(:,:), eps
  REAL(8) :: Fab(3), dist2, denominador, f_prod_q(SIZE(m))
  INTEGER :: a, b

  f_prod_q = 0.0
  DO a = 2, SIZE(m)
    DO b = 1, a - 1
      Fab = R(b,:) - R(a,:)
      dist2 = DOT_PRODUCT(Fab, Fab)
      denominador = SQRT(dist2 + eps*eps)

      ! forca
      Fab = G * m(a) * m(b) * Fab / (denominador*denominador*denominador)

      ! virial
      f_prod_q(a) = f_prod_q(a) + DOT_PRODUCT(R(a,:), Fab)
      f_prod_q(b) = f_prod_q(b) - DOT_PRODUCT(R(b,:), Fab)
    END DO
  END DO

END FUNCTION virial_potencial_amortecido_separado

FUNCTION termo_virial (G, m, R, soft) RESULT(termos)
  REAL(8), INTENT(IN) :: G, m(:), R(:,:)
  REAL(8), INTENT(IN), OPTIONAL :: soft
  REAL(8) :: eps, termos(2), dist2, den, termo, Fab(3)
  INTEGER :: a, b
  eps = 0.0
  IF (PRESENT(soft)) eps = soft
  termos = 0.0

  DO a = 2, SIZE(m)
    DO b = 1, a - 1
      Fab = R(b,:) - R(a,:)
      dist2 = DOT_PRODUCT(Fab, Fab)
      den = SQRT(dist2 + eps*eps)
      termo = G * m(a) * m(b) / den
      termos(1) = termos(1) - dist2 * termo / (den * den)
      termos(2) = termos(2) - termo
    END DO
  END DO
END FUNCTION termo_virial

FUNCTION tensor_inercia_individual (m, R) RESULT(tensor)
  REAL(8), INTENT(IN) :: m, R(3)
  REAL(8) :: tensor(3,3)
  INTEGER :: a, b

  DO a = 1, 3
    DO b = 1, 3
      IF (a /= b) tensor(a,b) = - m * R(a) * R(b)
    END DO
  END DO

  tensor(1,1) = m * (R(2)**2 + R(3)**2)
  tensor(2,2) = m * (R(1)**2 + R(3)**2)
  tensor(3,3) = m * (R(1)**2 + R(2)**2)
END FUNCTION tensor_inercia_individual

FUNCTION tensor_inercia_geral (ms, qs) RESULT(tensor)
  REAL(8), INTENT(IN) :: ms(:), qs(:,:)
  REAL(8), DIMENSION(3,3) :: tensor
  INTEGER :: a

  tensor(:,:) = 0.0

  DO a = 1, SIZE(ms)
    tensor = tensor + tensor_inercia_individual(ms(a), qs(a,:))
  END DO

END FUNCTION tensor_inercia_geral

END MODULE mecanica