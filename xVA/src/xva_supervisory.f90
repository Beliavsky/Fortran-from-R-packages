module xva_supervisory
  use trading, only : dp, str_len
  use xva_types, only : supervisory_cva_data_t, labeled_matrix_t, &
    relationship_correlation_t, commodity_risk_weight_t, equity_risk_weight_t
  use xva_math, only : same_text
  implicit none
  private

  public :: default_supervisory_cva_data
  public :: load_supervisory_cva_data
  public :: rating_weight
  public :: sector_risk_weight

contains

  subroutine default_supervisory_cva_data(data)
    type(supervisory_cva_data_t), intent(out) :: data
    integer :: i

    allocate(data%ir_eligible_tenors(5), data%ir_other_tenors(5), &
      data%ir_eligible_correlation(5,5), data%ir_other_correlation(5,5), &
      data%ir_risk_weight_eligible(5), data%ir_risk_weight_other(5))
    data%ir_eligible_tenors = [1.0_dp, 2.0_dp, 5.0_dp, 10.0_dp, 30.0_dp]
    data%ir_other_tenors = data%ir_eligible_tenors
    data%ir_eligible_correlation = reshape([ &
      1.00_dp, 0.91_dp, 0.72_dp, 0.55_dp, 0.31_dp, &
      0.91_dp, 1.00_dp, 0.87_dp, 0.72_dp, 0.45_dp, &
      0.72_dp, 0.87_dp, 1.00_dp, 0.91_dp, 0.68_dp, &
      0.55_dp, 0.72_dp, 0.91_dp, 1.00_dp, 0.83_dp, &
      0.31_dp, 0.45_dp, 0.68_dp, 0.83_dp, 1.00_dp ], [5,5])
    data%ir_other_correlation = 0.50_dp
    do i = 1, 5
      data%ir_other_correlation(i, i) = 1.0_dp
    end do
    data%ir_risk_weight_eligible = [0.0111_dp, 0.0093_dp, 0.0074_dp, &
      0.0074_dp, 0.0074_dp]
    data%ir_risk_weight_other = 0.0158_dp

    allocate(data%cs_tenors(5), data%cs_tenor_correlation(5,5))
    data%cs_tenors = [0.5_dp, 1.0_dp, 3.0_dp, 5.0_dp, 10.0_dp]
    data%cs_tenor_correlation = 0.90_dp
    do i = 1, 5
      data%cs_tenor_correlation(i, i) = 1.0_dp
    end do

    allocate(data%sectors(9), data%sector_risk_weight_ig(9), &
      data%sector_risk_weight_hy_nr(9))
    data%sectors = [character(len=str_len) :: &
      "Sovereigns including central banks and multilateral development banks", &
      "Local government, government-backed non-financials, education and public administration", &
      "Financials including government-backed financials", &
      "Basic materials, energy, industrials, agriculture, manufacturing, mining and quarrying", &
      "Consumer goods and services, transportation and storage, administrative and support service activities", &
      "Technology, telecommunications", &
      "Health care, utilities, professional and technical activities", &
      "Other sector", "Qualified Indices" ]
    data%sector_risk_weight_ig = [0.005_dp, 0.01_dp, 0.05_dp, 0.03_dp, &
      0.03_dp, 0.02_dp, 0.015_dp, 0.05_dp, 0.015_dp]
    data%sector_risk_weight_hy_nr = [0.02_dp, 0.04_dp, 0.12_dp, 0.07_dp, &
      0.085_dp, 0.055_dp, 0.05_dp, 0.12_dp, 0.05_dp]

    allocate(data%ratings(7), data%rating_weights(7))
    data%ratings = [character(len=str_len) :: "AAA", "AA", "A", "BBB", &
      "BB", "B", "CCC"]
    data%rating_weights = [0.007_dp, 0.007_dp, 0.008_dp, 0.010_dp, &
      0.020_dp, 0.030_dp, 0.100_dp]
  end subroutine default_supervisory_cva_data

  subroutine load_supervisory_cva_data(data_directory, data)
    character(len=*), intent(in) :: data_directory
    type(supervisory_cva_data_t), intent(out) :: data
    character(len=:), allocatable :: directory

    directory = trim(data_directory)
    if (len(directory) == 0) error stop "load_supervisory_cva_data: empty directory"
    if (directory(len(directory):len(directory)) /= "/" .and. &
        iachar(directory(len(directory):len(directory))) /= 92) then
      directory = directory // "/"
    end if

    call default_supervisory_cva_data(data)
    call load_correlation_matrix(directory // "IR_cormat_eligible_ccies.csv", &
      data%ir_eligible_tenors, data%ir_eligible_correlation)
    call load_correlation_matrix(directory // "IR_cormat_other_ccies.csv", &
      data%ir_other_tenors, data%ir_other_correlation)
    call load_correlation_matrix(directory // "CS_cormat_by_tenor.csv", &
      data%cs_tenors, data%cs_tenor_correlation)
    call load_ir_risk_weights(directory // "IR_RW.csv", &
      data%ir_risk_weight_eligible, data%ir_risk_weight_other)
    call load_sector_risk_weights(directory // "superv_risk_weights.csv", data)
    call load_rating_weights(directory // "RatingsMapping.csv", data)

    call load_labeled_matrix(directory // "IR_cormat_eligible_ccies.csv", &
      data%ir_eligible_full)
    call load_labeled_matrix(directory // "IR_cormat_other_ccies.csv", &
      data%ir_other_full)
    call load_labeled_matrix(directory // "CS_cormat_by_sector.csv", &
      data%cs_correlation_by_sector)
    call load_labeled_matrix(directory // "CS_cormat_by_tenor.csv", &
      data%cs_correlation_by_tenor_full)
    call load_labeled_matrix(directory // "CS_ref_Sector_cormat.csv", &
      data%cs_reference_sector_correlation)
    call load_labeled_matrix(directory // "CS_Sector_cpty_cormat.csv", &
      data%cs_sector_counterparty_correlation)
    call load_relationship_correlations(directory // "hedge_cpty_corr.csv", &
      data%hedge_counterparty_correlations)
    call load_commodity_risk_weights(directory // "COMM_RW.csv", &
      data%commodity_risk_weights)
    call load_equity_risk_weights(directory // "EQ_RW.csv", &
      data%equity_risk_weights)
  end subroutine load_supervisory_cva_data

  real(dp) function rating_weight(data, rating, found) result(value)
    type(supervisory_cva_data_t), intent(in) :: data
    character(len=*), intent(in) :: rating
    logical, intent(out), optional :: found
    integer :: i

    value = 0.0_dp
    if (present(found)) found = .false.
    do i = 1, size(data%ratings)
      if (same_text(data%ratings(i), rating)) then
        value = data%rating_weights(i)
        if (present(found)) found = .true.
        return
      end if
    end do
  end function rating_weight

  real(dp) function sector_risk_weight(data, sector, investment_grade, found) &
      result(value)
    type(supervisory_cva_data_t), intent(in) :: data
    character(len=*), intent(in) :: sector
    logical, intent(in) :: investment_grade
    logical, intent(out), optional :: found
    integer :: i

    value = 0.0_dp
    if (present(found)) found = .false.
    do i = 1, size(data%sectors)
      if (same_text(data%sectors(i), sector)) then
        if (investment_grade) then
          value = data%sector_risk_weight_ig(i)
        else
          value = data%sector_risk_weight_hy_nr(i)
        end if
        if (present(found)) found = .true.
        return
      end if
    end do
  end function sector_risk_weight


  subroutine load_labeled_matrix(filename, table)
    character(len=*), intent(in) :: filename
    type(labeled_matrix_t), intent(out) :: table
    character(len=8192) :: line
    character(len=str_len) :: fields(64)
    integer :: column_count
    integer :: field_count
    integer :: ios
    integer :: j
    integer :: row
    integer :: row_count
    integer :: unit

    open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "load_labeled_matrix: unable to open file"
    read(unit, '(a)', iostat=ios) line
    if (ios /= 0) error stop "load_labeled_matrix: unable to read header"
    call split_csv(line, fields, field_count)
    column_count = field_count - 1
    if (column_count <= 0) error stop "load_labeled_matrix: no numeric columns"

    row_count = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) > 0) row_count = row_count + 1
    end do
    rewind(unit)
    read(unit, '(a)', iostat=ios) line
    call split_csv(line, fields, field_count)

    allocate(table%row_labels(row_count), table%column_labels(column_count), &
      table%values(row_count, column_count))
    table%column_labels = fields(2:column_count + 1)
    row = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      call split_csv(line, fields, field_count)
      if (field_count < column_count + 1) cycle
      row = row + 1
      table%row_labels(row) = trim(fields(1))
      do j = 1, column_count
        table%values(row, j) = parse_percentage(fields(j + 1))
      end do
    end do
    close(unit)
    if (row /= row_count) error stop "load_labeled_matrix: unexpected row count"
  end subroutine load_labeled_matrix

  subroutine load_relationship_correlations(filename, values)
    character(len=*), intent(in) :: filename
    type(relationship_correlation_t), allocatable, intent(out) :: values(:)
    character(len=2048) :: line
    character(len=str_len) :: fields(8)
    integer :: count
    integer :: field_count
    integer :: ios
    integer :: row
    integer :: unit

    call count_data_rows(filename, count)
    allocate(values(count))
    open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "load_relationship_correlations: unable to open file"
    read(unit, '(a)', iostat=ios) line
    row = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      call split_csv(line, fields, field_count)
      if (field_count < 2) cycle
      row = row + 1
      values(row)%relationship = trim(fields(1))
      values(row)%correlation = parse_percentage(fields(2))
    end do
    close(unit)
    if (row /= count) error stop "load_relationship_correlations: unexpected row count"
  end subroutine load_relationship_correlations

  subroutine load_commodity_risk_weights(filename, values)
    character(len=*), intent(in) :: filename
    type(commodity_risk_weight_t), allocatable, intent(out) :: values(:)
    character(len=2048) :: line
    character(len=str_len) :: fields(8)
    integer :: count
    integer :: field_count
    integer :: ios
    integer :: row
    integer :: unit

    call count_data_rows(filename, count)
    allocate(values(count))
    open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "load_commodity_risk_weights: unable to open file"
    read(unit, '(a)', iostat=ios) line
    row = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      call split_csv(line, fields, field_count)
      if (field_count < 3) cycle
      row = row + 1
      values(row)%bucket_number = parse_integer(fields(1))
      values(row)%group = trim(fields(2))
      values(row)%risk_weight = parse_percentage(fields(3))
    end do
    close(unit)
    if (row /= count) error stop "load_commodity_risk_weights: unexpected row count"
  end subroutine load_commodity_risk_weights

  subroutine load_equity_risk_weights(filename, values)
    character(len=*), intent(in) :: filename
    type(equity_risk_weight_t), allocatable, intent(out) :: values(:)
    character(len=4096) :: line
    character(len=str_len) :: fields(12)
    integer :: count
    integer :: field_count
    integer :: ios
    integer :: row
    integer :: unit

    call count_data_rows(filename, count)
    allocate(values(count))
    open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "load_equity_risk_weights: unable to open file"
    read(unit, '(a)', iostat=ios) line
    row = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      call split_csv(line, fields, field_count)
      if (field_count < 5) cycle
      row = row + 1
      values(row)%bucket_number = parse_integer(fields(1))
      values(row)%size_class = trim(fields(2))
      values(row)%region = trim(fields(3))
      values(row)%sector = trim(fields(4))
      values(row)%risk_weight = parse_percentage(fields(5))
    end do
    close(unit)
    if (row /= count) error stop "load_equity_risk_weights: unexpected row count"
  end subroutine load_equity_risk_weights

  subroutine count_data_rows(filename, count)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: count
    character(len=8192) :: line
    integer :: ios
    integer :: unit

    open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "count_data_rows: unable to open file"
    read(unit, '(a)', iostat=ios) line
    count = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) > 0) count = count + 1
    end do
    close(unit)
  end subroutine count_data_rows

  subroutine load_correlation_matrix(filename, tenors, correlation)
    character(len=*), intent(in) :: filename
    real(dp), allocatable, intent(inout) :: tenors(:)
    real(dp), allocatable, intent(inout) :: correlation(:,:)
    character(len=4096) :: line
    character(len=str_len) :: fields(32)
    integer :: field_count
    integer :: ios
    integer :: row
    integer :: unit
    integer :: j

    open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "load_correlation_matrix: unable to open file"
    read(unit, '(a)', iostat=ios) line
    row = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      call split_csv(line, fields, field_count)
      if (field_count < size(tenors) + 1) cycle
      if (.not. is_numeric(fields(1))) cycle
      row = row + 1
      if (row > size(tenors)) exit
      tenors(row) = parse_number(fields(1))
      do j = 1, size(tenors)
        correlation(row, j) = parse_percentage(fields(j + 1))
      end do
    end do
    close(unit)
    if (row /= size(tenors)) error stop "load_correlation_matrix: unexpected row count"
  end subroutine load_correlation_matrix

  subroutine load_ir_risk_weights(filename, eligible, other)
    character(len=*), intent(in) :: filename
    real(dp), allocatable, intent(inout) :: eligible(:)
    real(dp), allocatable, intent(inout) :: other(:)
    character(len=1024) :: line
    character(len=str_len) :: fields(8)
    integer :: field_count
    integer :: ios
    integer :: row
    integer :: unit

    open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "load_ir_risk_weights: unable to open file"
    read(unit, '(a)', iostat=ios) line
    row = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      call split_csv(line, fields, field_count)
      if (field_count < 3 .or. .not. is_numeric(fields(1))) cycle
      row = row + 1
      if (row > size(eligible)) exit
      eligible(row) = parse_percentage(fields(2))
      other(row) = parse_percentage(fields(3))
    end do
    close(unit)
    if (row /= size(eligible)) error stop "load_ir_risk_weights: unexpected row count"
  end subroutine load_ir_risk_weights

  subroutine load_sector_risk_weights(filename, data)
    character(len=*), intent(in) :: filename
    type(supervisory_cva_data_t), intent(inout) :: data
    character(len=4096) :: line
    character(len=str_len) :: fields(8)
    integer :: field_count
    integer :: ios
    integer :: row
    integer :: unit

    open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "load_sector_risk_weights: unable to open file"
    read(unit, '(a)', iostat=ios) line
    row = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      call split_csv(line, fields, field_count)
      if (field_count < 3) cycle
      row = row + 1
      if (row > size(data%sectors)) exit
      data%sectors(row) = trim(fields(1))
      data%sector_risk_weight_ig(row) = parse_percentage(fields(2))
      data%sector_risk_weight_hy_nr(row) = parse_percentage(fields(3))
    end do
    close(unit)
    if (row /= size(data%sectors)) then
      error stop "load_sector_risk_weights: unexpected row count"
    end if
  end subroutine load_sector_risk_weights

  subroutine load_rating_weights(filename, data)
    character(len=*), intent(in) :: filename
    type(supervisory_cva_data_t), intent(inout) :: data
    character(len=1024) :: line
    character(len=str_len) :: fields(4)
    integer :: field_count
    integer :: ios
    integer :: row
    integer :: unit

    open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "load_rating_weights: unable to open file"
    read(unit, '(a)', iostat=ios) line
    row = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      call split_csv(line, fields, field_count)
      if (field_count < 2) cycle
      row = row + 1
      if (row > size(data%ratings)) exit
      data%ratings(row) = trim(adjustl(fields(1)))
      data%rating_weights(row) = parse_percentage(fields(2))
    end do
    close(unit)
    if (row /= size(data%ratings)) error stop "load_rating_weights: unexpected row count"
  end subroutine load_rating_weights

  subroutine split_csv(line, fields, field_count)
    character(len=*), intent(in) :: line
    character(len=str_len), intent(out) :: fields(:)
    integer, intent(out) :: field_count
    character(len=str_len) :: current
    logical :: quoted
    integer :: i
    integer :: position

    fields = ""
    current = ""
    field_count = 1
    position = 0
    quoted = .false.
    i = 1
    do while (i <= len_trim(line))
      select case (line(i:i))
      case ('"')
        if (quoted .and. i < len_trim(line) .and. line(i + 1:i + 1) == '"') then
          position = position + 1
          if (position <= len(current)) current(position:position) = '"'
          i = i + 1
        else
          quoted = .not. quoted
        end if
      case (',')
        if (quoted) then
          position = position + 1
          if (position <= len(current)) current(position:position) = ','
        else
          if (field_count <= size(fields)) fields(field_count) = trim(current)
          field_count = field_count + 1
          current = ""
          position = 0
        end if
      case default
        if (iachar(line(i:i)) /= 13) then
          position = position + 1
          if (position <= len(current)) current(position:position) = line(i:i)
        end if
      end select
      i = i + 1
    end do
    if (field_count <= size(fields)) fields(field_count) = trim(current)
    field_count = min(field_count, size(fields))
  end subroutine split_csv

  logical function is_numeric(text) result(value)
    character(len=*), intent(in) :: text
    real(dp) :: parsed
    integer :: ios

    read(text, *, iostat=ios) parsed
    value = ios == 0
  end function is_numeric

  real(dp) function parse_percentage(text) result(value)
    character(len=*), intent(in) :: text
    character(len=str_len) :: work
    integer :: index_value

    work = trim(adjustl(text))
    index_value = index(work, "%")
    if (index_value > 0) then
      work(index_value:index_value) = " "
      value = parse_number(work) / 100.0_dp
    else
      value = parse_number(work)
    end if
  end function parse_percentage

  integer function parse_integer(text) result(value)
    character(len=*), intent(in) :: text
    integer :: ios

    read(text, *, iostat=ios) value
    if (ios /= 0) error stop "parse_integer: invalid integer field"
  end function parse_integer

  real(dp) function parse_number(text) result(value)
    character(len=*), intent(in) :: text
    integer :: ios

    read(text, *, iostat=ios) value
    if (ios /= 0) error stop "parse_number: invalid numeric field"
  end function parse_number

end module xva_supervisory
