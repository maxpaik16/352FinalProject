
# File for building neural network to model stability

# necessary imports
import csv
import torch
from torch import nn
from torch.utils.data import DataLoader, Dataset


# define class for planetary data
# each entry includes the positions and velocity of each planet and the
# stability of the system after 1000 Neptune periods
class PlanetaryDataset(Dataset):
    def __init__(self, csvfile:
        # define constructor
        self.data = []
        with open(csvfile) as datafile:
            csvReader = csv.reader(datafile)
            i = 0
            for row in csvReader:
                if int(row[-1]) == 1 or i % 4 == 0:
                    # add to balance stable and unstable data
                    self.data.append(torch.FloatTensor(list(map(lambda x : float(x), row))))
                i += 1

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        # get data at given inset
        planets = self.data[idx][:-1]
        return {"data": planets, "stable":torch.FloatTensor([self.data[idx][-1]])}


# define neural network class
class StabilityNetwork(nn.Module):
    def __init__(self):
        super(StabilityNetwork, self).__init__()

        def init_weights(m):
            if type(m) == nn.Linear:
                torch.nn.init.xavier_uniform_(m.weight)
                m.bias.data.fill_(0.01)

        # choose shallow network for relatively simple task
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(8*4, 128),
            nn.ReLU(),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, 1)
        )
        self.linear_relu_stack = self.linear_relu_stack.apply(init_weights)

    def forward(self, x):
        logits = self.linear_relu_stack(x)
        return logits

# function for training the model on a given set of data
def train_loop(dataloader, model, loss_fn, optimizer):
    size = len(dataloader.dataset)
    for batch, xy in enumerate(dataloader):
        for i in range((len(xy['data']))):
            # Compute prediction and loss
            pred = model(xy['data'][i])
            loss = loss_fn(pred, xy['stable'][i])

            # Backpropagation
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            if batch % 100 == 0:
                loss, current = loss.item(), i * len(xy['data'][i])
                print(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")

# function for evaluating the model's performance on a given set of data
def test_loop(dataloader, model, loss_fn):
    size = len(dataloader.dataset)
    test_loss, correct = 0, 0

    with torch.no_grad():
        for batch, xy in enumerate(dataloader):
            for i in range((len(xy['data']))):
                pred = model(xy['data'][i])
                test_loss += loss_fn(pred, xy['stable'][i]).item()

                # check binary classification
                res = int(pred[0].item() > 0)
                correct += int(res == int(xy['stable'][i].item()))

    test_loss /= size
    correct /= size
    print(f"Test Error: \n Accuracy: {(100*correct):>0.1f}%, Avg loss: {test_loss:>8f} \n")
    return test_loss


if __name__ == '__main__':
    # build model
    device = 'cpu'
    model = StabilityNetwork().to(device)

    # set up training data and hyperparameters
    training_data = PlanetaryDataset('training.csv')
    batch_size = 16
    training_loader = DataLoader(training_data, batch_size=batch_size, shuffle=True)
    learning_rate = .1
    epochs = 5
    loss_fn = nn.BCEWithLogitsLoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, betas=(0.9, 0.999))

    # run main training/evaluation loop
    losses = []
    for t in range(epochs):
        print(f"Epoch {t+1}\n-------------------------------")
        model.train()
        train_loop(training_loader, model, loss_fn, optimizer)
        model.eval()
        l = test_loop(training_loader, model, loss_fn)
        losses.append(l)
    print("Done!")
