#include "refcom.hpp"

using namespace CommandLineProcessing;

int parse_args(int argc, char* argv[], InputArgs& in_args){

  ArgvParser cmd;
  const char* op_arg = "operation";
  const char* rf_arg = "ref_file_name";
  const char* f1_arg = "rd1_file_name";
  const char* f2_arg = "rd2_file_name";
  const char* f3_arg = "com_file_name";
  const char* tc_arg = "thread_count";
  const char* l1_arg = "rd1_length";
  const char* l2_arg = "rd2_length";

  cmd.setIntroductoryDescription("Reference-based Compression and Decompression");
  cmd.addErrorCode(0, "Success");
  cmd.addErrorCode(1, "Error");

  cmd.setHelpOption("h", "help", "Print this help page");

  cmd.defineOption(op_arg, "Operation : compression/decompression",
                   ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
  //cmd.defineOptionAlternative(op_arg, "");

  cmd.defineOption(rf_arg, "Reference file name",
                   ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
  //cmd.defineOptionAlternative(rf_arg, "");

  cmd.defineOption(f1_arg, "Read-1 file name",
                   ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
  //cmd.defineOptionAlternative(f1_arg, "");

  cmd.defineOption(f2_arg, "Read-2 file name",
                   ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
  //cmd.defineOptionAlternative(f2_arg, "");

  cmd.defineOption(f3_arg, "Compressed data file name",
                   ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
  //cmd.defineOptionAlternative(f3_arg, "");

  cmd.defineOption(tc_arg, "Thread count",
                   ArgvParser::OptionRequiresValue);
  //cmd.defineOptionAlternative(tc_arg, "");

  cmd.defineOption(l1_arg, "Read-1 length",
                   ArgvParser::OptionRequiresValue);
  //cmd.defineOptionAlternative(l1_arg, "");

  cmd.defineOption(l2_arg, "Read-2 length",
                   ArgvParser::OptionRequiresValue);
  //cmd.defineOptionAlternative(l2_arg, "");

  int result = cmd.parse(argc, argv);

  //Make sure we get the right command line args
  if (result != ArgvParser::NoParserError) {
    std::cerr << cmd.parseErrorDescription(result) << std::endl;
    return 1;
  }

  if(cmd.foundOption(op_arg)) {
      in_args.operation = cmd.optionValue(op_arg);
  } else {
      std::cerr << "Required option missing: " << op_arg << std::endl;
      return 1;
  }

  if(cmd.foundOption(rf_arg)) {
      in_args.refFileName = cmd.optionValue(rf_arg);
  } else {
      std::cerr << "Required option missing: " << rf_arg << std::endl;
      return 1;
  }

  if(cmd.foundOption(f1_arg)) {
      in_args.rd1FileName = cmd.optionValue(f1_arg);
  } else {
      std::cerr << "Required option missing: " << f1_arg << std::endl;
      return 1;
  }

  if(cmd.foundOption(f2_arg)) {
      in_args.rd2FileName = cmd.optionValue(f2_arg);
  } else {
      std::cerr << "Required option missing: " << f2_arg << std::endl;
      return 1;
  }

  if(cmd.foundOption(f3_arg)) {
      in_args.comFileName = cmd.optionValue(f3_arg);
  } else {
      std::cerr << "Required option missing: " << f3_arg << std::endl;
      return 1;
  }

  if(cmd.foundOption(tc_arg)) {
      in_args.threadCount = (std::uint32_t) std::stoi(cmd.optionValue(tc_arg));
  } else {
      std::cerr << "Required option missing: " << tc_arg << std::endl;
      return 1;
  }

  if(cmd.foundOption(l1_arg)) {
      in_args.rd1Length = (std::uint32_t) std::stoi(cmd.optionValue(l1_arg));
  } else {
      std::cerr << "Required option missing: " << l1_arg << std::endl;
      return 1;
  }

  if(cmd.foundOption(l2_arg)) {
      in_args.rd2Length = (std::uint32_t) std::stoi(cmd.optionValue(l2_arg));
  } else {
      std::cerr << "Required option missing: " << l2_arg << std::endl;
      return 1;
  }

  std::cerr << "-----------------------------------------------------" << std::endl;
  std::cerr << "Operation          : " << in_args.operation   << std::endl;
  std::cerr << "Ref File Name      : " << in_args.refFileName << std::endl;
  std::cerr << "Rd1 File Name      : " << in_args.rd1FileName << std::endl;
  std::cerr << "Rd2 File Name      : " << in_args.rd2FileName << std::endl;
  std::cerr << "Com File Name      : " << in_args.comFileName << std::endl;
  std::cerr << "Thread Count       : " << in_args.threadCount << std::endl;
  std::cerr << "Rd1 Length         : " << in_args.rd1Length   << std::endl;
  std::cerr << "Rd2 Length         : " << in_args.rd2Length   << std::endl;
  std::cerr << "-----------------------------------------------------" << std::endl;

  return 0;
}

int refcom_compression(InputArgs& in_args, CompressionDataStructures& comDS) {
    if (in_args.threadCount < 2) {
        std::cerr << "Thread count must be at least 2" << std::endl;
        return 1;
    }
    if ((in_args.rd1Length < 1 ) || (in_args.rd1Length > 158)) {
        std::cerr << "Rd1 length must be between 1 and 158" << std::endl;
        return 1;
    }
    if ((in_args.rd2Length < 1 ) || (in_args.rd2Length > 158)) {
        std::cerr << "Rd2 length must be between 1 and 158" << std::endl;
        return 1;
    }
    
    build_index(in_args, comDS);
    
    align_reads(in_args, comDS);
    
    //compress_reads(in_args, comDS);
    //TODO: Move below inside compress_reads
    free(comDS.ref_bases); comDS.ref_length = 0; //Free memory
    
    return 0;
}

int refcom_decompression(InputArgs& in_args, DecompressionDataStructures& decomDS) {
    
    return 0;
}

int main(int argc, char *argv[]) {

    //Parse command line arguments
    InputArgs cargs;
    
    if(parse_args(argc, argv, cargs))
        return 1;
    
    CompressionDataStructures comDS;
    DecompressionDataStructures decomDS;
    
    if (cargs.operation.compare("compression")) {
        if (refcom_compression(cargs, comDS))
            return 1;
    } else if (cargs.operation.compare("decompression")) {
        if (refcom_decompression(cargs, decomDS))
            return 1;
    } else {
        std::cerr << "Unrecognized operation : " << cargs.operation << std::endl;
        std::cerr << "Must be compression or decompression" << std::endl;
        return 1;
    }
    
    return 0;
}























































