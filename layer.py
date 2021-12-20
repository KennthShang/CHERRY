import  torch
from    torch import nn
from    torch.nn import functional as F
from    utils import sparse_dropout, dot


class GraphConvolution(nn.Module):


    def __init__(self, node_dim,input_dim, output_dim, num_features_nonzero,
                 dropout=0.,
                 is_sparse_inputs=False,
                 bias=False,
                 activation = F.relu,
                 featureless=False,
                 learn_weight=True):
        super(GraphConvolution, self).__init__()


        self.dropout = dropout
        self.bias = bias
        self.activation = activation
        self.is_sparse_inputs = is_sparse_inputs
        self.featureless = featureless
        self.num_features_nonzero = num_features_nonzero

        self.weight = nn.Parameter(torch.randn(input_dim, output_dim))
        self.edge_w = nn.Parameter(torch.rand(node_dim, node_dim))
        self.bias = None
        self.learn_weight=True
        if bias:
            self.bias = nn.Parameter(torch.zeros(output_dim))


    def forward(self, inputs):
        # print('inputs:', inputs)
        x, support = inputs

        if self.training:
            x = F.dropout(x, self.dropout)

        xw = torch.mm(x, self.weight)
        support_w = support
        #support_w = torch.mm(support, self.edge_w)
        #support_w = torch.mul(support, self.edge_w)
        #if self.learn_weight:
            #support_w = torch.mm(support, self.edge_w)
        #    support_w = torch.mul(support, self.edge_w)
        #else:
        #    support_w = support
        out = torch.mm(support_w , xw)

        if self.bias is not None:
            out += self.bias

        return self.activation(out), support

