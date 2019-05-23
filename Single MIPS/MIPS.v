// Single Cycle MIPS
//=========================================================
// Input/Output Signals:
// positive-edge triggered         clk
// active low asynchronous reset   rst_n
// instruction memory interface    IR_addr, IR
// output for testing purposes     RF_writedata  
//=========================================================
// Wire/Reg Specifications:
// control signals             MemToReg, MemRead, MemWrite, 
//                             RegDST, RegWrite, Branch, 
//                             Jump, ALUSrc, ALUOp
// ALU control signals         ALUctrl
// ALU input signals           ALUin1, ALUin2
// ALU output signals          ALUresult, ALUzero
// instruction specifications  r, j, jal, jr, lw, sw, beq
// sign-extended signal        SignExtend
// MUX output signals          MUX_RegDST, MUX_MemToReg, 
//                             MUX_Src, MUX_Branch, MUX_Jump
// registers input signals     Reg_R1, Reg_R2, Reg_W, WriteData 
// registers                   Register
// registers output signals    ReadData1, ReadData2
// data memory contral signals CEN, OEN, WEN
// data memory output signals  ReadDataMem
// program counter/address     PCin, PCnext, JumpAddr, BranchAddr
//=========================================================

module SingleCycle_MIPS( 
    clk,
    rst_n,
    IR_addr,
    IR,
    RF_writedata,
    ReadDataMem,
    CEN,
    WEN,
    A,
    ReadData2,
    OEN
);

//==== in/out declaration =================================
    //-------- processor ----------------------------------
    input         clk, rst_n;
    input  [31:0] IR;
    output [31:0] IR_addr, RF_writedata;
    //-------- data memory --------------------------------
    input  [31:0] ReadDataMem;  // read_data from memory
    output        CEN;  // chip_enable, 0 when you read/write data from/to memory
    output        WEN;  // write_enable, 0 when you write data into SRAM & 1 when you read data from SRAM
    output  [6:0] A;  // address
    output [31:0] ReadData2;  // write_data to memory
    output        OEN;  // output_enable, 0

//==== reg/wire declaration ===============================

    reg [31:0] IR_addr, RF_writedata;
// control signals 
    wire MemRead, MemWrite, RegWrite, Branch, Jump, ALUSrc;
    wire [1:0] ALUOp, MemToReg, RegDST;

// ALU control signals  
    wire [2:0] ALUctrl;
    wire JRctrl;

// ALU input signals 
    reg [31:0] ALUin1, ALUin2;

// ALU output signals 
    reg [31:0] ALUresult, ALUzero;

// sign-extended signal
    wire [31:0] SignExtend;

// MUX output signals 
    wire [31:0] MUX_MemToReg, MUX_Src, MUX_Branch, MUX_Jump, MUX_Jr;
    wire [4:0] MUX_RegDST;

// registers input signals
    reg [4:0] Reg_R1, Reg_R2, Reg_W;

// registers
    reg [31:0] Register [31:0];

// registers output signals
    reg [31:0] ReadData1, ReadData2, WriteData;

// data memory contral signals
    wire WEN, CEN, OEN;

// data memory output signals
    wire [31:0] ReadDataMem;

// program counter/address
    wire [31:0] PCin, PCnext, JumpAddr, BranchAddr;

    wire BranchDeside;
    wire [31:0] BranchShift2;


//==== combinational part =================================


// control unit
    assign RegDST[0] = (~IR[31]) & (~IR[30]) & (~IR[29]) & (~IR[28]) & (~IR[27]) & (~IR[26]);
    assign Jump = ((~IR[31]) & (~IR[30]) & (~IR[29]) & (~IR[28]) & (IR[27]) & (~IR[26])) | ((~IR[31]) & (~IR[30]) & (~IR[29]) & (~IR[28]) & (IR[27]) & (IR[26]));
	assign Branch = (~IR[31]) & (~IR[30]) & (~IR[29]) & (IR[28]) & (~IR[27]) & (~IR[26]);
	assign MemRead = (IR[31]) & (~IR[30]) & (~IR[29]) & (~IR[28]) & (IR[27]) & (IR[26]);
	assign MemToReg[0] = (IR[31]) & (~IR[30]) & (~IR[29]) & (~IR[28]) & (IR[27]) & (IR[26]);
	assign MemWrite = (IR[31]) & (~IR[30]) & (IR[29]) & (~IR[28]) & (IR[27]) & (IR[26]);
	assign ALUSrc  = ((IR[31]) & (~IR[30]) & (~IR[29]) & (~IR[28]) & (IR[27]) & (IR[26])) | ((IR[31]) & (~IR[30]) & (IR[29]) & (~IR[28]) & (IR[27]) & (IR[26])); 
	assign RegWrite = ((~IR[31]) & (~IR[30]) & (~IR[29]) & (~IR[28]) & (~IR[27]) & (~IR[26])) | ((IR[31]) & (~IR[30]) & (~IR[29]) & (~IR[28]) & (IR[27]) & (IR[26])) | ((~IR[31]) & (~IR[30]) & (~IR[29]) & (~IR[28]) & (IR[27]) & (IR[26]));
	assign ALUOp[1] = ((~IR[31]) & (~IR[30]) & (~IR[29]) & (~IR[28]) & (~IR[27]) & (~IR[26]));
	assign ALUOp[0] = ((~IR[31]) & (~IR[30]) & (~IR[29]) & (IR[28]) & (~IR[27]) & (~IR[26]));
    assign RegDST[1] = (~IR[31]) & (~IR[30]) & (~IR[29]) & (~IR[28]) & (IR[27]) & (IR[26]);
    assign MemToReg[1] = (~IR[31]) & (~IR[30]) & (~IR[29]) & (~IR[28]) & (IR[27]) & (IR[26]);
    assign CEN = (~MemWrite) & (~MemRead);
    assign WEN = (~MemWrite) & (MemRead);
    assign OEN = 0;

// alu
	assign ALUctrl[2] = ((ALUOp[0]) | (ALUOp[1] & IR[1]));
	assign ALUctrl[1] = ((~ALUOp[1]) | (~IR[2]));
	assign ALUctrl[0] = ((ALUOp[1]) & (IR[0] | IR[3]));
	assign JRctrl = (ALUOp[1] & ~IR[5] & ~IR[4] & IR[3] & ~IR[2] & ~IR[1] & ~IR[0]);
    

    always@(*)
    begin
        ALUin1 = ReadData1;
        ALUin2 = MUX_Src;

        case(ALUctrl)
            3'b000: begin
                ALUzero = 0;
                ALUresult = ALUin1 & ALUin2;
            end
            3'b001: begin
                ALUzero = 0;
                ALUresult = ALUin1 | ALUin2;
            end
            3'b010: begin
                ALUzero = 0;
                ALUresult = ALUin1 + ALUin2;
            end
            3'b110: begin
                ALUresult = ALUin1 - ALUin2;
                if (ALUresult == 0) 
                    begin
                        ALUzero = 1;
                    end
            end
            3'b111: begin
                ALUzero = 0;
                ALUresult = (ALUin1 - ALUin2) >> 31;
            end
            default: begin
                ALUzero = 0;
                ALUresult = ALUin1;
            end
        endcase
    end

    assign A = ALUresult[8:2];

// pc
    assign PCnext = IR_addr + 32'b0100;
    assign SignExtend = {{14{IR[15]}}, IR[15:0]};
    assign BranchAddr = {SignExtend[29:0], 2'b00};
    assign JumpAddr = {PCnext[31:28], {IR[25:0], 2'b00}};
    //assign MUX_RegDST = RegDST ? IR[15:11] : IR[20:16];
    assign MUX_RegDST = RegDST[1] ? {5'b11111} : (RegDST ? IR[15:11] : IR[20:16]);
    assign MUX_Src = ALUSrc ? SignExtend : ReadData2;
    //assign MUX_MemToReg = MemToReg ? ReadDataMem : ALUresult;
    assign MUX_MemToReg = MemToReg[1] ? IR_addr + 31'b100 : (MemToReg[0] ? ReadDataMem : ALUresult);
    assign MUX_Branch = BranchDeside ? PCnext + BranchAddr : PCnext;
    assign MUX_Jump = Jump ? JumpAddr : MUX_Branch;
    assign MUX_Jr = JRctrl ? ReadData1 : MUX_Jump;
    assign PCin = MUX_Jr;
    and (BranchDeside, ALUzero, Branch);


// register file
    always@(*)
    begin
        Reg_R1 = IR[25:21];
        Reg_R2 = IR[20:16];
        ReadData1 = Register[Reg_R1];
        ReadData2 = Register[Reg_R2];
        Reg_W = MUX_RegDST;
        WriteData = MUX_MemToReg;
        RF_writedata = WriteData;
    end



//==== sequential part ====================================

// pc
	always@(posedge clk or negedge rst_n)
    begin
        if (rst_n == 1'b0)
        begin
        	IR_addr <= 0;
        end
        else
        begin
			IR_addr <= PCin;
        end
    end
// register file
	integer i;
	always@(posedge clk or negedge rst_n)
    begin
        if (rst_n == 1'b0) 
        begin
            for (i = 0; i < 32; i = i + 1)
            begin
                Register[i] <= 32'b0;
            end
        end
        else if (RegWrite == 1'b1) Register[MUX_RegDST] <= WriteData;
    end

//=========================================================
endmodule
