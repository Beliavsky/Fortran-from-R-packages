program convert_prices
   use comparison_asset_data, only : asset_price_data, read_asset_prices, &
      write_asset_prices_binary, asset_data_ok
   implicit none

   type(asset_price_data) :: prices
   integer :: status
   character(len=512) :: csv_file, binary_file
   character(len=256) :: message

   call get_command_argument(1, csv_file)
   call get_command_argument(2, binary_file)
   if (len_trim(csv_file) == 0) csv_file = '../../asset_class_etf_prices.csv'
   if (len_trim(binary_file) == 0) binary_file = '../../asset_class_etf_prices.bin'

   call read_asset_prices(trim(csv_file), prices, status, message)
   if (status /= asset_data_ok) error stop trim(message)
   call write_asset_prices_binary(trim(binary_file), prices, status, message)
   if (status /= asset_data_ok) error stop trim(message)
   write(*, '(a,i0,a,i0,a)') 'Wrote ', size(prices%prices, 1), ' rows and ', &
      size(prices%prices, 2), ' assets to '//trim(binary_file)
end program convert_prices
